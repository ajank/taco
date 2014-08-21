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

StatsSet::StatsSet(PositionWeightMatrix *M0, PositionWeightMatrix *M1, int margin)
{
  StatsSet::M0 = M0;
  StatsSet::M1 = M1;
  stats.min_offset = stats.min_offset_same = stats.min_offset_opposite = -M1->length - margin;
  stats.max_offset = stats.max_offset_same = stats.max_offset_opposite = M0->length + margin;
  offset_shift_same = offset_shift_opposite = 0;
  initialize_carriers();
}

StatsSet::StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin)
{
  StatsSet::M0 = h1->parts[0].M;
  StatsSet::M1 = h2->parts[0].M;

  stats.min_offset_same = h1->pair_start - h2->pair_end - margin;
  stats.min_offset_opposite = h1->pair_start - (h2->parts[0].M->length - h2->pair_start) - margin;
  stats.min_offset = min(stats.min_offset_same, stats.min_offset_opposite);

  stats.max_offset_same = h1->pair_end - h2->pair_start + margin;
  stats.max_offset_opposite = h1->pair_end - (h2->parts[0].M->length - h2->pair_end) + margin;
  stats.max_offset = max(stats.max_offset_same, stats.max_offset_opposite);

  offset_shift_same = h2->pair_start - h1->pair_start;
  offset_shift_opposite = (h2->parts[0].M->length - h2->pair_end) - h1->pair_start;

  initialize_carriers();
}

void StatsSet::initialize_carriers()
{
  stats.same_carrier.resize(stats.max_offset - stats.min_offset + 1);
  stats.same = &stats.same_carrier[0] - stats.min_offset;
  stats.opposite_carrier.resize(stats.max_offset - stats.min_offset + 1);
  stats.opposite = &stats.opposite_carrier[0] - stats.min_offset;
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
  vector<PositionWeightMatrix::MotifMatch>::iterator M0_it = M0->matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::iterator M0_matches_end = M0->matches.end();
  vector<PositionWeightMatrix::MotifMatch>::iterator M1_it = M1->matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::iterator M1_matches_end = M1->matches.end();

  for (list<NarrowPeak::Region>::const_iterator jt = regions->regions.begin(); jt != regions->regions.end(); jt++)
  {
    // taking pseudo_start and pseudo_end, the frequencies will be comparable regardles of the region fragmentation
    int jt_pseudo_start = jt->start - M0->length / 2;
    int jt_pseudo_end = jt->end - M0->length / 2;
    while (M0_it != M0_matches_end && M0_it->chrom_id < jt->chrom_id)
      M0_it++;
    while (M0_it != M0_matches_end && M0_it->chrom_id == jt->chrom_id && M0_it->start < jt_pseudo_start)
      M0_it++;
    while (M0_it != M0_matches_end && M0_it->chrom_id == jt->chrom_id && M0_it->start < jt_pseudo_end)
      M0_filtered_matches.push_back(*M0_it++);

    jt_pseudo_start = jt->start - M1->length / 2;
    jt_pseudo_end = jt->end - M1->length / 2;
    while (M1_it != M1_matches_end && M1_it->chrom_id < jt->chrom_id)
      M1_it++;
    while (M1_it != M1_matches_end && M1_it->chrom_id == jt->chrom_id && M1_it->start < jt_pseudo_start)
      M1_it++;
    while (M1_it != M1_matches_end && M1_it->chrom_id == jt->chrom_id && M1_it->start < jt_pseudo_end)
      M1_filtered_matches.push_back(*M1_it++);
  }
}

void StatsSet::calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M0_matches, const vector<PositionWeightMatrix::MotifMatch> &M1_matches)
{
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
  }
}

void StatsSet::returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M0_paired_matches, int offset, bool same_orientation) const
{
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
