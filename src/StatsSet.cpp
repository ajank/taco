/*
    StatsSet.cpp

    Copyright (C) 2011-2013  Aleksander Jankowski <ajank@mimuw.edu.pl>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <stdlib.h>

#include "StatsSet.h"

void Stats::initialize_hits_carriers()
{
  hits_same_carrier.resize(max_offset - min_offset + 1);
  hits_same = &hits_same_carrier[0] - min_offset;
  hits_opposite_carrier.resize(max_offset - min_offset + 1);
  hits_opposite = &hits_opposite_carrier[0] - min_offset;
}

void Stats::initialize_stats_carriers()
{
  same_carrier.resize(max_offset - min_offset + 1);
  same = &same_carrier[0] - min_offset;
  opposite_carrier.resize(max_offset - min_offset + 1);
  opposite = &opposite_carrier[0] - min_offset;
}

StatsSet::StatsSet(const vector<PositionWeightMatrix *> Mi, int margin)
{
  arity = Mi.size();
  if (arity < 2)
  {
    cerr << "StatsSet: arity < 2" << endl;
    exit(1);
  }

  StatsSet::Mi = Mi;
  Mi_filtered_matches.resize(arity);
  min_offset.resize(arity - 1);
  max_offset.resize(arity - 1);
  initialize_stats(stats, 1, 0, Mi[0]->length, margin);
  for (int level = 1; level < arity; level++)
  {
    int i = level - 1;
    cerr << "min_offset " << min_offset[i] << " max_offset " << max_offset[i] << endl;
  }
}

StatsSet::StatsSet(PositionWeightMatrix *M0, PositionWeightMatrix *M1, int margin)
{
  arity = 2;
  Mi.push_back(M0);
  Mi.push_back(M1);
  Mi_filtered_matches.resize(arity);

  stats.min_offset = stats.min_offset_same = stats.min_offset_opposite = -M1->length - margin;
  min_offset.push_back(stats.min_offset);

  stats.max_offset = stats.max_offset_same = stats.max_offset_opposite = M0->length + margin;
  max_offset.push_back(stats.max_offset);

  offset_shift_same = offset_shift_opposite = 0;

  stats.initialize_hits_carriers();
}

StatsSet::StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin)
{
  arity = 2;
  Mi.push_back(h1->parts[0].M);
  Mi.push_back(h2->parts[0].M);
  Mi_filtered_matches.resize(arity);

  stats.min_offset_same = h1->pair_start - h2->pair_end - margin;
  stats.min_offset_opposite = h1->pair_start - (h2->parts[0].M->length - h2->pair_start) - margin;
  stats.min_offset = min(stats.min_offset_same, stats.min_offset_opposite);
  min_offset.push_back(stats.min_offset);

  stats.max_offset_same = h1->pair_end - h2->pair_start + margin;
  stats.max_offset_opposite = h1->pair_end - (h2->parts[0].M->length - h2->pair_end) + margin;
  stats.max_offset = max(stats.max_offset_same, stats.max_offset_opposite);
  max_offset.push_back(stats.max_offset);

  offset_shift_same = h2->pair_start - h1->pair_start;
  offset_shift_opposite = (h2->parts[0].M->length - h2->pair_end) - h1->pair_start;

  stats.initialize_hits_carriers();
}

void StatsSet::initialize_stats(Stats &stats, int level, int pair_start, int pair_end, int margin)
{
  PositionWeightMatrix *M1 = Mi[level];

  stats.min_offset = stats.min_offset_same = stats.min_offset_opposite = pair_start - M1->length - margin;
  min_offset[level - 1] = min(min_offset[level - 1], stats.min_offset);

  stats.max_offset = stats.max_offset_same = stats.max_offset_opposite = pair_end + margin;
  max_offset[level - 1] = max(max_offset[level - 1], stats.max_offset);

  offset_shift_same = offset_shift_opposite = 0;

  if (level + 1 == arity)
    stats.initialize_hits_carriers();
  else
  {
    stats.initialize_stats_carriers();

    for (int offset = stats.min_offset_same; offset <= stats.max_offset_same; offset++)
      initialize_stats(stats.same[offset], level + 1, min(pair_start, offset), max(pair_end, M1->length + offset), margin);

    for (int offset = stats.min_offset_opposite; offset <= stats.max_offset_opposite; offset++)
      initialize_stats(stats.opposite[offset], level + 1, min(pair_start, offset), max(pair_end, M1->length + offset), margin);
  }
}

void StatsSet::addHit(const vector<int> &offset, const vector<bool> &same_orientation)
{
  if ((int) offset.size() + 1 != arity)
  {
    cerr << "StatsSet: addHit: arity " << arity << ", but got offset length " << offset.size() << endl;
    exit(1);
  }

  if ((int) same_orientation.size() + 1 != arity)
  {
    cerr << "StatsSet: addHit: arity " << arity << ", but got same_orientation length " << same_orientation.size() << endl;
    exit(1);
  }

  Stats s = stats;
  for (int level = 1; level < arity; level++)
  {
    int i = level - 1;
    if (s.min_offset <= offset[i] && offset[i] <= s.max_offset)
    {
      if (level + 1 == arity)
      {
        if (same_orientation[i])
          s.hits_same[offset[i]]++;
        else
          s.hits_opposite[offset[i]]++;
      }
      else
      {
        if (same_orientation[i])
          s = s.same[offset[i]];
        else
          s = s.opposite[offset[i]];
      }
    }
    else
      break;
  }
}

void StatsSet::addHit(vector<int> offset, vector<bool> same_orientation, int level, const vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> &jrt_head, const vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> &jrt_tail)
{
/*  cerr << "addHit";
  for (int i = 0; i < offset.size(); i++) cerr << " " << offset[i];
  cerr << " / ";
  for (int i = 0; i < same_orientation.size(); i++) cerr << " " << same_orientation[i];
  cerr << " / " << level << " / ";
  //for (int i = 0; i < jrt_head.size(); i++) cerr << " " << jrt_head[i];
  cerr << " / ";
  //for (int i = 0; i < jrt_tail.size(); i++) cerr << " " << jrt_tail[i];
  cerr << endl;*/

  if ((int) same_orientation.size() + 1 != level)
  {
    cerr << "StatsSet: same_orientation length " << same_orientation.size() << " + 1 != level " << level << endl;
    exit(1);
  }

  if ((int) offset.size() + 1 != level)
  {
    cerr << "StatsSet: offset length " << offset.size() << " + 1 != level " << level << endl;
    exit(1);
  }

  if ((int) jrt_head.size() != arity)
  {
    cerr << "StatsSet: jrt_head length " << jrt_head.size() << " != arity " << arity << endl;
    exit(1);
  }

  if ((int) jrt_tail.size() != arity)
  {
    cerr << "StatsSet: jrt_tail length " << jrt_tail.size() << " != arity " << arity << endl;
    exit(1);
  }

  offset.push_back(int());
  same_orientation.push_back(bool());

  for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head[level]; jrt != jrt_tail[level]; jrt++)
  {
    if (jrt_head[0]->strand == '+')
    {
      offset.back() = jrt->start - jrt_head[0]->start;
      same_orientation.back() = jrt->strand == '+';
    }
    else
    {
      offset.back() = jrt_head[0]->start + Mi[0]->length - jrt->start - Mi[level]->length;
      same_orientation.back() = jrt->strand != '+';
    }

    if (level + 1 == arity)
      addHit(offset, same_orientation);
    else
      addHit(offset, same_orientation, level + 1, jrt_head, jrt_tail);
  }
}

void StatsSet::addHit(int offset, bool same_orientation)
{
  if (stats.min_offset <= offset && offset <= stats.max_offset)
  {
    if (same_orientation)
      stats.hits_same[offset]++;
    else
      stats.hits_opposite[offset]++;
  }
}

long int StatsSet::getStats(int offset, bool same_orientation) const
{
  if (stats.min_offset <= offset && offset <= stats.max_offset)
  {
    if (same_orientation)
      return stats.hits_same[offset];
    else
      return stats.hits_opposite[offset];
  }
  else
  {
    cerr << "StatsSet: offset " << offset << " beyond range [" << stats.min_offset << ", " << stats.max_offset << "]" << endl;
    exit(1);
  }
}

long int StatsSet::getMaximumStats(int *returned_offset, bool *returned_same_orientation) const
{
  long int hits = -1;

  for (int offset = stats.min_offset_same; offset <= stats.max_offset_same; offset++)
    if (stats.hits_same[offset] > hits)
    {
      hits = stats.hits_same[offset];
      *returned_offset = offset + offset_shift_same;
      *returned_same_orientation = true;
    }

  for (int offset = stats.min_offset_opposite; offset <= stats.max_offset_opposite; offset++)
    if (stats.hits_opposite[offset] > hits)
    {
      hits = stats.hits_opposite[offset];
      *returned_offset = offset + offset_shift_opposite;
      *returned_same_orientation = false;
    }

  return hits;
}

void StatsSet::addMotifMatchesFromRegions(const NarrowPeak *regions)
{
  vector<vector<PositionWeightMatrix::MotifMatch>::iterator> Mi_it, Mi_matches_end;
  for (int i = 0; i < arity; i++)
  {
    Mi_it.push_back(Mi[i]->matches.begin());
    Mi_matches_end.push_back(Mi[i]->matches.end());
  }

  for (list<NarrowPeak::Region>::const_iterator jt = regions->regions.begin(); jt != regions->regions.end(); jt++)
  {
    for (int i = 0; i < arity; i++)
    {
      // taking pseudo_start and pseudo_end, the frequencies will be comparable regardles of the region fragmentation
      int jt_pseudo_start = jt->start - Mi[i]->length / 2;
      int jt_pseudo_end = jt->end - Mi[i]->length / 2;
      while (Mi_it[i] != Mi_matches_end[i] && Mi_it[i]->chrom_id < jt->chrom_id)
        Mi_it[i]++;
      while (Mi_it[i] != Mi_matches_end[i] && Mi_it[i]->chrom_id == jt->chrom_id && Mi_it[i]->start < jt_pseudo_start)
        Mi_it[i]++;
      while (Mi_it[i] != Mi_matches_end[i] && Mi_it[i]->chrom_id == jt->chrom_id && Mi_it[i]->start < jt_pseudo_end)
        Mi_filtered_matches[i].push_back(*Mi_it[i]++);
    }
  }
}

void StatsSet::calculateStats(const vector<vector<PositionWeightMatrix::MotifMatch> > &Mi_matches)
{
  if ((int) Mi_matches.size() != arity)
  {
    cerr << "StatsSet: calculateStats: arity " << arity << ", but got Mi_matches length " << Mi_matches.size() << endl;
    exit(1);
  }

  vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> jrt_head, jrt_tail;
  vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> Mi_matches_end;
  for (int i = 0; i < arity; i++)
  {
    jrt_head.push_back(Mi_matches[i].begin());
    jrt_tail.push_back(Mi_matches[i].begin());
    Mi_matches_end.push_back(Mi_matches[i].end());
  }

/* see the inequalities for dimer case below */

  for (; jrt_head[0] != Mi_matches_end[0]; jrt_head[0]++)
  {
    for (int i = 1; i < arity; i++)
    {
      while (jrt_head[i] != Mi_matches_end[i] && jrt_head[i]->chrom_id < jrt_head[0]->chrom_id)
        jrt_head[i]++;
      while (jrt_head[i] != Mi_matches_end[i] && jrt_head[i]->chrom_id == jrt_head[0]->chrom_id && jrt_head[i]->start < jrt_head[0]->start + min_offset[i - 1])
        jrt_head[i]++;

      while (jrt_tail[i] != Mi_matches_end[i] && jrt_tail[i]->chrom_id < jrt_head[0]->chrom_id)
        jrt_tail[i]++;
      while (jrt_tail[i] != Mi_matches_end[i] && jrt_tail[i]->chrom_id == jrt_head[0]->chrom_id && jrt_tail[i]->start <= jrt_head[0]->start + max_offset[i - 1])
        jrt_tail[i]++;
    }

    addHit(vector<int>(), vector<bool>(), 1, jrt_head, jrt_tail);
  }
}

void StatsSet::calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M0_matches, const vector<PositionWeightMatrix::MotifMatch> &M1_matches)
{
  PositionWeightMatrix *M0 = Mi[0];
  PositionWeightMatrix *M1 = Mi[1];

  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_head = M1_matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_tail = M1_matches.begin();

/* We need to filter jrts. There are two cases:
(1) same orientation: looking for jrts satisfying
      stats.min_offset <= jrt->start - irt->start <= stats.max_offset
      irt->start + stats.min_offset <= jrt->start <= irt->start + stats.max_offset

(2) opposite orientation: looking for jrts satisfying
      stats.min_offset <= irt->end - jrt->end <= stats.max_offset
      stats.min_offset <= irt->start + M0->length - jrt->start - M1->length <= stats.max_offset
      - stats.max_offset <= - irt->start - M0->length + jrt->start + M1->length <= - stats.min_offset
      irt->start + M0->length - M1->length - stats.max_offset <= jrt->start <= irt->start + M0->length - M1->length - stats.min_offset
      irt->start + M0->length - M1->length - stats.max_offset <= jrt->start <= irt->start + M0->length - M1->length - stats.min_offset
    but
      M0->length - M1->length - stats.max_offset = M0->length - M1->length - M0->length - margin + 1 = stats.min_offset
      M0->length - M1->length - stats.min_offset = M0->length - M1->length + M1->length + margin - 1 = stats.max_offset
    so the above inequalities are equivalent to
      irt->start + stats.min_offset <= jrt->start <= irt->start + stats.max_offset. */

  for (vector<PositionWeightMatrix::MotifMatch>::const_iterator irt = M0_matches.begin(); irt != M0_matches.end(); irt++)
  {
    while (jrt_head != M1_matches.end() && jrt_head->chrom_id < irt->chrom_id)
      jrt_head++;
    while (jrt_head != M1_matches.end() && jrt_head->chrom_id == irt->chrom_id && jrt_head->start < irt->start + stats.min_offset)
      jrt_head++;

    while (jrt_tail != M1_matches.end() && jrt_tail->chrom_id < irt->chrom_id)
      jrt_tail++;
    while (jrt_tail != M1_matches.end() && jrt_tail->chrom_id == irt->chrom_id && jrt_tail->start <= irt->start + stats.max_offset)
      jrt_tail++;

    for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
    {
      if (irt->strand == '+')
        addHit(jrt->start - irt->start, jrt->strand == '+');
      else
        addHit(irt->start + M0->length - jrt->start - M1->length, jrt->strand != '+');
    }

/*    vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> new_jrt_head;
    new_jrt_head.push_back(irt);
    new_jrt_head.push_back(jrt_head);
    vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> new_jrt_tail;
    new_jrt_tail.push_back(irt);
    new_jrt_tail.push_back(jrt_tail);
    addHit(vector<int>(), vector<bool>(), 1, new_jrt_head, new_jrt_tail);*/
  }
}

void StatsSet::returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M0_paired_matches, int offset, bool same_orientation) const
{
  //FIXME
  PositionWeightMatrix *M0 = Mi[0];
  PositionWeightMatrix *M1 = Mi[1];
  vector<PositionWeightMatrix::MotifMatch> M0_filtered_matches = Mi_filtered_matches[0];
  vector<PositionWeightMatrix::MotifMatch> M1_filtered_matches = Mi_filtered_matches[1];

  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_head = M1_filtered_matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_tail = M1_filtered_matches.begin();

  for (vector<PositionWeightMatrix::MotifMatch>::const_iterator irt = M0_filtered_matches.begin(); irt != M0_filtered_matches.end(); irt++)
  {
    while (jrt_head != M1_filtered_matches.end() && jrt_head->chrom_id < irt->chrom_id)
      jrt_head++;
    while (jrt_head != M1_filtered_matches.end() && jrt_head->chrom_id == irt->chrom_id && jrt_head->start < irt->start + stats.min_offset)
      jrt_head++;

    while (jrt_tail != M1_filtered_matches.end() && jrt_tail->chrom_id < irt->chrom_id)
      jrt_tail++;
    while (jrt_tail != M1_filtered_matches.end() && jrt_tail->chrom_id == irt->chrom_id && jrt_tail->start <= irt->start + stats.max_offset)
      jrt_tail++;

    if (irt->strand == '+')
    {
      for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
        if (jrt->start - irt->start == offset && (jrt->strand == '+') == same_orientation)
          M0_paired_matches.push_back(*irt);
    }
    else
    {
      for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
        if (irt->start + M0->length - jrt->start - M1->length == offset && (jrt->strand != '+') == same_orientation)
          M0_paired_matches.push_back(*irt);
    }
  }
}
