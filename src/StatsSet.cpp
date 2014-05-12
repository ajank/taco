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

StatsSet::StatsSet(PositionWeightMatrix *M1, PositionWeightMatrix *M2, int margin)
{
  StatsSet::M1 = M1;
  StatsSet::M2 = M2;
  min_offset = min_offset_same = min_offset_opposite = -M2->length - margin;
  max_offset = max_offset_same = max_offset_opposite = M1->length + margin;
  offset_shift_same = offset_shift_opposite = 0;
  initialize_carriers();
}

StatsSet::StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin)
{
  StatsSet::M1 = h1->M1;
  StatsSet::M2 = h2->M1;

  min_offset_same = h1->pair_start - h2->pair_end - margin;
  min_offset_opposite = h1->pair_start - (h2->M1->length - h2->pair_start) - margin;
  min_offset = min(min_offset_same, min_offset_opposite);

  max_offset_same = h1->pair_end - h2->pair_start + margin;
  max_offset_opposite = h1->pair_end - (h2->M1->length - h2->pair_end) + margin;
  max_offset = max(max_offset_same, max_offset_opposite);

  offset_shift_same = h2->pair_start - h1->pair_start;
  offset_shift_opposite = (h2->M1->length - h2->pair_end) - h1->pair_start;

  initialize_carriers();
}

StatsSet::StatsSet(const StatsSet &other)
{
  copy_constants(other);
  initialize_carriers();
  copy(other.same_carrier, other.same_carrier + max_offset - min_offset + 1, same_carrier);
  copy(other.opposite_carrier, other.opposite_carrier + max_offset - min_offset + 1, opposite_carrier);
}

StatsSet &StatsSet::operator= (const StatsSet &other)
{
  copy_constants(other);

  Stats *same_carrier_temp = new Stats[max_offset - min_offset + 1];
  copy(other.same_carrier, other.same_carrier + max_offset - min_offset + 1, same_carrier_temp);
  delete[] same_carrier;
  same_carrier = same_carrier_temp;
  same = same_carrier - min_offset;

  Stats *opposite_carrier_temp = new Stats[max_offset - min_offset + 1];
  copy(other.opposite_carrier, other.opposite_carrier + max_offset - min_offset + 1, opposite_carrier_temp);
  delete[] opposite_carrier;
  opposite_carrier = opposite_carrier_temp;
  opposite = opposite_carrier - min_offset;

  return *this;
}

void StatsSet::initialize_carriers()
{
  same_carrier = new Stats[max_offset - min_offset + 1];
  same = same_carrier - min_offset;
  opposite_carrier = new Stats[max_offset - min_offset + 1];
  opposite = opposite_carrier - min_offset;
}

void StatsSet::copy_constants(const StatsSet &other)
{
  M1 = other.M1;
  M2 = other.M2;
  min_offset = other.min_offset;
  min_offset_same = other.min_offset_same;
  min_offset_opposite = other.min_offset_opposite;
  max_offset = other.max_offset;
  max_offset_same = other.max_offset_same;
  max_offset_opposite = other.max_offset_opposite;
  offset_shift_same = other.offset_shift_same;
  offset_shift_opposite = other.offset_shift_opposite;
}

void StatsSet::addHit(int offset, bool same_orientation)
{
  if (min_offset <= offset && offset <= max_offset)
  {
    if (same_orientation)
      same[offset].addHit();
    else
      opposite[offset].addHit();
  }
}

Stats *StatsSet::getStats(int offset, bool same_orientation) const
{
  if (min_offset <= offset && offset <= max_offset)
  {
    if (same_orientation)
      return &same[offset];
    else
      return &opposite[offset];
  }
  else
  {
    cerr << "StatsSet: offset " << offset << " beyond range [" << min_offset << ", " << max_offset << "]" << endl;
    exit(1);
  }
}

long int StatsSet::getMaximumStats(int *returned_offset, bool *returned_same_orientation) const
{
  long int hits = -1;

  for (int offset = min_offset_same; offset <= max_offset_same; offset++)
    if (same[offset].getHits() > hits)
    {
      hits = same[offset].getHits();
      *returned_offset = offset + offset_shift_same;
      *returned_same_orientation = true;
    }

  for (int offset = min_offset_opposite; offset <= max_offset_opposite; offset++)
    if (opposite[offset].getHits() > hits)
    {
      hits = opposite[offset].getHits();
      *returned_offset = offset + offset_shift_opposite;
      *returned_same_orientation = false;
    }

  return hits;
}

void StatsSet::addMotifMatchesFromRegions(const NarrowPeak *regions)
{
  vector<PositionWeightMatrix::MotifMatch>::iterator M1_it = M1->matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::iterator M1_matches_end = M1->matches.end();
  vector<PositionWeightMatrix::MotifMatch>::iterator M2_it = M2->matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::iterator M2_matches_end = M2->matches.end();

  for (list<NarrowPeak::Region>::const_iterator jt = regions->regions.begin(); jt != regions->regions.end(); jt++)
  {
    // taking pseudo_start and pseudo_end, the frequencies will be comparable regardles of the region fragmentation
    int jt_pseudo_start = jt->start - M1->length / 2;
    int jt_pseudo_end = jt->end - M1->length / 2;
    while (M1_it != M1_matches_end && M1_it->chrom_id < jt->chrom_id)
      M1_it++;
    while (M1_it != M1_matches_end && M1_it->chrom_id == jt->chrom_id && M1_it->start < jt_pseudo_start)
      M1_it++;
    while (M1_it != M1_matches_end && M1_it->chrom_id == jt->chrom_id && M1_it->start < jt_pseudo_end)
      M1_filtered_matches.push_back(*M1_it++);

    jt_pseudo_start = jt->start - M2->length / 2;
    jt_pseudo_end = jt->end - M2->length / 2;
    while (M2_it != M2_matches_end && M2_it->chrom_id < jt->chrom_id)
      M2_it++;
    while (M2_it != M2_matches_end && M2_it->chrom_id == jt->chrom_id && M2_it->start < jt_pseudo_start)
      M2_it++;
    while (M2_it != M2_matches_end && M2_it->chrom_id == jt->chrom_id && M2_it->start < jt_pseudo_end)
      M2_filtered_matches.push_back(*M2_it++);
  }
}

void StatsSet::calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M1_matches, const vector<PositionWeightMatrix::MotifMatch> &M2_matches)
{
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_head = M2_matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_tail = M2_matches.begin();

/* We need to filter jrts. There are two cases:
(1) same orientation: looking for jrts satisfying
      min_offset <= jrt->start - irt->start <= max_offset
      irt->start + min_offset <= jrt->start <= irt->start + max_offset

(2) opposite orientation: looking for jrts satisfying
      min_offset <= irt->end - jrt->end <= max_offset
      min_offset <= irt->start + M1->length - jrt->start - M2->length <= max_offset
      - max_offset <= - irt->start - M1->length + jrt->start + M2->length <= - min_offset
      irt->start + M1->length - M2->length - max_offset <= jrt->start <= irt->start + M1->length - M2->length - min_offset
      irt->start + M1->length - M2->length - max_offset <= jrt->start <= irt->start + M1->length - M2->length - min_offset
    but
      M1->length - M2->length - max_offset = M1->length - M2->length - M1->length - margin + 1 = min_offset
      M1->length - M2->length - min_offset = M1->length - M2->length + M2->length + margin - 1 = max_offset
    so the above inequalities are equivalent to
      irt->start + min_offset <= jrt->start <= irt->start + max_offset. */

  for (vector<PositionWeightMatrix::MotifMatch>::const_iterator irt = M1_matches.begin(); irt != M1_matches.end(); irt++)
  {
    while (jrt_head != M2_matches.end() && jrt_head->chrom_id < irt->chrom_id)
      jrt_head++;
    while (jrt_head != M2_matches.end() && jrt_head->chrom_id == irt->chrom_id && jrt_head->start < irt->start + min_offset)
      jrt_head++;

    while (jrt_tail != M2_matches.end() && jrt_tail->chrom_id < irt->chrom_id)
      jrt_tail++;
    while (jrt_tail != M2_matches.end() && jrt_tail->chrom_id == irt->chrom_id && jrt_tail->start <= irt->start + max_offset)
      jrt_tail++;

    for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
    {
      if (irt->strand == '+')
        addHit(jrt->start - irt->start, jrt->strand == '+');
      else
        addHit(irt->start + M1->length - jrt->start - M2->length, jrt->strand != '+');
    }
  }
}

void StatsSet::returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M1_paired_matches, int offset, bool same_orientation) const
{
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_head = M2_filtered_matches.begin();
  vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt_tail = M2_filtered_matches.begin();

  for (vector<PositionWeightMatrix::MotifMatch>::const_iterator irt = M1_filtered_matches.begin(); irt != M1_filtered_matches.end(); irt++)
  {
    while (jrt_head != M2_filtered_matches.end() && jrt_head->chrom_id < irt->chrom_id)
      jrt_head++;
    while (jrt_head != M2_filtered_matches.end() && jrt_head->chrom_id == irt->chrom_id && jrt_head->start < irt->start + min_offset)
      jrt_head++;

    while (jrt_tail != M2_filtered_matches.end() && jrt_tail->chrom_id < irt->chrom_id)
      jrt_tail++;
    while (jrt_tail != M2_filtered_matches.end() && jrt_tail->chrom_id == irt->chrom_id && jrt_tail->start <= irt->start + max_offset)
      jrt_tail++;

    if (irt->strand == '+')
    {
      for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
        if (jrt->start - irt->start == offset && (jrt->strand == '+') == same_orientation)
          M1_paired_matches.push_back(*irt);
    }
    else
    {
      for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jrt = jrt_head; jrt != jrt_tail; jrt++)
        if (irt->start + M1->length - jrt->start - M2->length == offset && (jrt->strand != '+') == same_orientation)
          M1_paired_matches.push_back(*irt);
    }
  }
}

StatsSet::~StatsSet()
{
  delete[] same_carrier;
  delete[] opposite_carrier;
}
