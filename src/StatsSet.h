/*
    StatsSet.h

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

#ifndef _STATS_SET_H
#define _STATS_SET_H

#include "HypothesesSet.h"
#include "NarrowPeak.h"
#include "PositionWeightMatrix.h"

using namespace std;

class Stats
{
  public:
    // first of all, know the limits:
    int min_offset, min_offset_same, min_offset_opposite;
    int max_offset, max_offset_same, max_offset_opposite;

    // then: either count the number of hits if level == arity...
    long int *hits_same, *hits_opposite;
    vector<long int> hits_same_carrier, hits_opposite_carrier;
    void initialize_hits_carriers();

    // ...or go one level down:
    Stats *same, *opposite;
    vector<Stats> same_carrier, opposite_carrier;
    void initialize_stats_carriers();
};

class StatsSet
{
  public:
    StatsSet(const vector<PositionWeightMatrix *> Mi, int margin);
    StatsSet(PositionWeightMatrix *M0, PositionWeightMatrix *M1, int margin);
    StatsSet(const vector<HypothesesSet::Hypothesis *> hi, int margin);
    StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin);

    int arity; // 2 for dimers, 3 for trimers etc.
    vector<PositionWeightMatrix *> Mi;
    vector<vector<PositionWeightMatrix::MotifMatch> > Mi_filtered_matches;
    int offset_shift_same, offset_shift_opposite;
    vector<int> min_offset, max_offset;
    Stats stats;

    void addHit(const vector<int> &offset, const vector<bool> &same_orientation);
    void addHit(vector<int> offset, vector<bool> same_orientation, int level, const vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> &jrt_head, const vector<vector<PositionWeightMatrix::MotifMatch>::const_iterator> &jrt_tail);
    void addHit(int offset, bool same_orientation);
    long int getStats(int offset, bool same_orientation) const;
    long int getMaximumStats(int *returned_offset, bool *returned_same_orientation) const;
    void calculateStats(const vector<vector<PositionWeightMatrix::MotifMatch> > &Mi_matches);
    void calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M0_matches, const vector<PositionWeightMatrix::MotifMatch> &M1_matches);

    void addMotifMatchesFromRegions(const NarrowPeak *regions);
    void calculateStats();
    void returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M0_paired_matches, int offset, bool same_orientation) const;

  private:
    void initialize_stats(Stats &stats, int level, int pair_start, int pair_end, int margin);
};

inline void StatsSet::calculateStats()
{
  return calculateStats(Mi_filtered_matches);
}

#endif
