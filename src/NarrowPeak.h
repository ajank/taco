/*
    NarrowPeak.h

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

#ifndef _NARROW_PEAK_H
#define _NARROW_PEAK_H

#include <list>
#include <map>
#include <string>

#include "Genome.h"

using namespace std;

#define RESTRICTION_NO_PEAK 1
#define RESTRICTION_NO_SIGNAL_VALUE 2

const int max_motif_midpoint_offset = 1000;

class NarrowPeak
{
  public:
    NarrowPeak();

    struct Region
    {
      int dataset_id, chrom_id;
      int start, end;
      double signalValue;
      int peak;
      bool operator<(const Region &other) const;
    };
    list<Region> regions;
    int restrictions;

    long int location_dist[max_motif_midpoint_offset + 1]; // for each (absolute value of) offset between motif midpoints, how many possible locations for such a motif complex are there?

    void readNarrowPeakFile(int dataset_id, const string &peakFile);
    void update_location_dist(const Genome &fa);
    void removeRegionsWithMaskedPeaks(const Genome &fa);
    void removeMostlyMaskedRegions(const Genome &fa);
    void restrictToTopSignalRegions(int limit, bool each_dataset_separately);
    void setRegionSize(int size);
    void mergeOverlappingRegions(bool each_dataset_separately);
  private:
    static bool sortFunc_signalValue(const Region i, const Region j);
    static bool sortFunc_dataset_id_signalValue(const Region i, const Region j);
};

#endif
