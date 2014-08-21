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
    long int *hits_same, *hits_opposite; // either count the number of hits if we reached the arity...
    vector<long int> hits_same_carrier, hits_opposite_carrier;

    Stats *same, *opposite; // ...or go one step down, increasing the arity
    vector<Stats> same_carrier, opposite_carrier;

    // in both of the cases, follow the limits:
    int min_offset, min_offset_same, min_offset_opposite;
    int max_offset, max_offset_same, max_offset_opposite;
};

class StatsSet
{
  public:
    StatsSet(PositionWeightMatrix *M0, PositionWeightMatrix *M1, int margin);
    StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin);

    PositionWeightMatrix *M0, *M1;
    vector<PositionWeightMatrix::MotifMatch> M0_filtered_matches, M1_filtered_matches;
    int offset_shift_same, offset_shift_opposite;
    Stats stats;

    void addHit(int offset, bool same_orientation);
    long int getStats(int offset, bool same_orientation) const;
    long int getMaximumStats(int *returned_offset, bool *returned_same_orientation) const;
    void calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M0_matches, const vector<PositionWeightMatrix::MotifMatch> &M1_matches);

    void addMotifMatchesFromRegions(const NarrowPeak *regions);
    void calculateStats();
    void returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M0_paired_matches, int offset, bool same_orientation) const;

  private:
    void initialize_carriers();
};

inline void StatsSet::calculateStats()
{
  return calculateStats(M0_filtered_matches, M1_filtered_matches);
}

#endif
