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
    Stats();
    void addHit();
    long int getHits() const;

  private:
    long int hits;
};

class StatsSet
{
  public:
    StatsSet(PositionWeightMatrix *M0, PositionWeightMatrix *M1, int margin);
    StatsSet(HypothesesSet::Hypothesis *h1, HypothesesSet::Hypothesis *h2, int margin);

    PositionWeightMatrix *M0, *M1;
    vector<PositionWeightMatrix::MotifMatch> M0_filtered_matches, M1_filtered_matches;
    int min_offset, min_offset_same, min_offset_opposite;
    int max_offset, max_offset_same, max_offset_opposite;
    int offset_shift_same, offset_shift_opposite;

    void addHit(int offset, bool same_orientation);
    Stats *getStats(int offset, bool same_orientation) const;
    long int getMaximumStats(int *returned_offset, bool *returned_same_orientation) const;
    void calculateStats(const vector<PositionWeightMatrix::MotifMatch> &M0_matches, const vector<PositionWeightMatrix::MotifMatch> &M1_matches);

    void addMotifMatchesFromRegions(const NarrowPeak *regions);
    void calculateStats();
    void returnPairedMatches(vector<PositionWeightMatrix::MotifMatch> &M0_paired_matches, int offset, bool same_orientation) const;

  private:
    void initialize_carriers();
    void copy_constants(const StatsSet &other);
    Stats *same, *opposite;
    vector<Stats> same_carrier, opposite_carrier;
};

inline Stats::Stats()
{
  hits = 0;
}

inline void Stats::addHit()
{
  hits++;
}

inline long int Stats::getHits() const
{
  return hits;
}

inline void StatsSet::calculateStats()
{
  return calculateStats(M0_filtered_matches, M1_filtered_matches);
}

#endif
