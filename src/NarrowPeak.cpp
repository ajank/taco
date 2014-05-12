/*
    NarrowPeak.cpp

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

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <zlib.h>

#include "NarrowPeak.h"

NarrowPeak::NarrowPeak()
{
  for (int i = 0; i < max_motif_midpoint_offset + 1; i++)
    location_dist[i] = -1;
}

void NarrowPeak::readNarrowPeakFile(int dataset_id, const string &peakFile)
{
  cout << "Reading narrowPeak/BED file: " << peakFile << endl;
  gzFile f = gzopen(peakFile.c_str(), "r");
  if (f == NULL)
  {
    cerr << peakFile << ": " << strerror(errno) << endl;
    exit(1);
  }

  struct Region row;
  char line[1024], chrom[1024];
  int lineNum = 0, numFields = 0;

  while (gzgets(f, line, 1024))
  {
    lineNum++;

    // guess the file format: ENCODE narrowPeak (BED 6+4, includes signalValue), BED 6 (includes score) or BED 3
    if (!numFields || numFields == 10)
    {
      if (sscanf(line, "%s%d%d %*s %*d %*c %lf%*f%*f%d", chrom, &row.start, &row.end, &row.signalValue, &row.peak) == 5)
      {
        numFields = 10;
        if (row.peak == -1)
        {
          row.peak = (row.end - row.start) / 2;
          restrictions |= RESTRICTION_NO_PEAK;
        }
      }
      else if (numFields)
      {
        cerr << peakFile << ": not enough values in line " << lineNum << " (expected narrowPeak format)" << endl;
        exit(1);
      }
    }

    if (!numFields || numFields == 6)
    {
      if (sscanf(line, "%s%d%d %*s %lf %*c", chrom, &row.start, &row.end, &row.signalValue) == 4) // take score as signalValue
      {
        numFields = 6;
        row.peak = (row.end - row.start) / 2;
        restrictions |= RESTRICTION_NO_PEAK;
      }
      else if (numFields)
      {
        cerr << peakFile << ": not enough values in line " << lineNum << " (expected BED 6 format)" << endl;
        exit(1);
      }
    }

    if (!numFields || numFields == 3)
    {
      if (sscanf(line, "%s%d%d", chrom, &row.start, &row.end) == 3)
      {
        numFields = 3;
        row.peak = (row.end - row.start) / 2;
        row.signalValue = -1;
        restrictions |= RESTRICTION_NO_PEAK | RESTRICTION_NO_SIGNAL_VALUE;
      }
      else if (numFields)
      {
        cerr << peakFile << ": not enough values in line " << lineNum << " (expected BED 3 format)" << endl;
        exit(1);
      }
    }

    if (!numFields)
    {
      cerr << peakFile << ": looks like neither narrowPeak nor BED format" << endl;
      exit(1);
    }

    row.chrom_id = Genome::get_chrom_id(chrom);
    row.dataset_id = dataset_id;
    regions.push_back(row);
  }

  gzclose(f);
}

bool NarrowPeak::Region::operator<(const Region &other) const
{
  if (this->chrom_id < other.chrom_id) return true;
  else if (this->chrom_id > other.chrom_id) return false;
  else return this->start < other.start;
}

void NarrowPeak::update_location_dist(const Genome &fa)
{
  for (int i = 0; i < max_motif_midpoint_offset + 1; i++)
    location_dist[i] = 0;

  int old_chrom_id = -1;
  map<int, vector<char> >::const_iterator jt;
  deque<pair<int, int> > sr; // subregions

  for (list<Region>::iterator it = regions.begin(); it != regions.end(); it++)
  {
    if (old_chrom_id != it->chrom_id)
    {
      old_chrom_id = it->chrom_id;
      jt = fa.chroms.find(it->chrom_id);
      if (jt == fa.chroms.end())
      {
        cerr << "NarrowPeak: chromosome \"" << Genome::chrom_names[it->chrom_id] << "\" not found in the genome" << endl;
        exit(1);
      }

      sr.clear();
    }

    int start = -1, end = it->start;
    while (end < it->end)
    {
      if (Genome::isUnmaskedOrSingletonN(jt->second, end))
      {
        if (start == -1) start = end; // start a new subregion

        if (!sr.empty()) // iterate over all partner midpoint locations in the previous subregions
          for (deque<pair<int, int> >::iterator si = sr.begin(); si != sr.end(); /*si++*/)
          {
            for (int i = min(end - si->first, max_motif_midpoint_offset); i > end - si->second; i--)
              location_dist[i]++;
            if (end - si->second > max_motif_midpoint_offset) si = sr.erase(si); else si++;
          }
      }
      else
      {
        if (start != -1)
        {
          for (int i = 0; i < min(end - start, max_motif_midpoint_offset + 1); i++)
            location_dist[i] += end - start - i;
          sr.push_back(make_pair(start, end));
        }
        start = -1;
      }

      end++;
    }

    if (start != -1)
    {
      for (int i = 0; i < min(end - start, max_motif_midpoint_offset + 1); i++)
        location_dist[i] += end - start - i;
      sr.push_back(make_pair(start, end));
    }
  }
}

void NarrowPeak::removeRegionsWithMaskedPeaks(const Genome &fa)
{
  if (restrictions & RESTRICTION_NO_PEAK)
    cout << "NarrowPeak: peak positions are missing for some regions, MaskedRegions=Peak must be applied carefully" << endl;

  int old_chrom_id = -1;
  map<int, vector<char> >::const_iterator jt;

  for (list<Region>::iterator it = regions.begin(); it != regions.end(); /*it++*/)
  {
    if (old_chrom_id != it->chrom_id)
    {
      old_chrom_id = it->chrom_id;
      jt = fa.chroms.find(it->chrom_id);
      if (jt == fa.chroms.end())
      {
        cerr << "NarrowPeak: chromosome \"" << Genome::chrom_names[it->chrom_id] << "\" not found in the genome" << endl;
        exit(1);
      }
    }

    if (Genome::isUnmaskedOrSingletonN(jt->second, it->start + it->peak))
      it++;
    else
      regions.erase(it++);
  }
}

void NarrowPeak::removeMostlyMaskedRegions(const Genome &fa)
{
  int old_chrom_id = -1;
  map<int, vector<char> >::const_iterator jt;

  for (list<Region>::iterator it = regions.begin(); it != regions.end(); /*it++*/)
  {
    if (old_chrom_id != it->chrom_id)
    {
      old_chrom_id = it->chrom_id;
      jt = fa.chroms.find(it->chrom_id);
      if (jt == fa.chroms.end())
      {
        cerr << "NarrowPeak: chromosome \"" << Genome::chrom_names[it->chrom_id] << "\" not found in the genome" << endl;
        exit(1);
      }
    }

    int count = 0;
    for (int i = it->start; i < it->end; i++)
      if (Genome::isUnmaskedOrSingletonN(jt->second, i)) count++;

    if (2 * count > it->end - it->start)
      it++;
    else
      regions.erase(it++);
  }
}

bool NarrowPeak::sortFunc_signalValue(const Region i, const Region j)
{
  return i.signalValue > j.signalValue;
}

bool NarrowPeak::sortFunc_dataset_id_signalValue(const Region i, const Region j)
{
  if (i.dataset_id < j.dataset_id) return true;
  else if (i.dataset_id > j.dataset_id) return false;
  else return i.signalValue > j.signalValue;
}

void NarrowPeak::restrictToTopSignalRegions(int limit, bool each_dataset_separately)
{
  if (restrictions & RESTRICTION_NO_SIGNAL_VALUE)
  {
    cerr << "NarrowPeak: signal values are missing for some regions, RegionCount=* cannot be applied" << endl;
    exit(1);
  }

  if (each_dataset_separately)
  {
    regions.sort(sortFunc_dataset_id_signalValue);

    int old_dataset_id = -1, count = 0;
    for (list<Region>::iterator it = regions.begin(); it != regions.end(); /*it++*/)
    {
      count++;
      if (old_dataset_id != it->dataset_id) count = 1;
      if (count == 1) old_dataset_id = it->dataset_id;

      if (count > limit) regions.erase(it++); else it++;
    }
  }
  else
  {
    regions.sort(sortFunc_signalValue);

    int count = 0;
    for (list<Region>::iterator it = regions.begin(); it != regions.end(); /*it++*/)
      if (++count > limit) regions.erase(it++); else it++;
  }
}

void NarrowPeak::setRegionSize(int size)
{
  for (list<Region>::iterator it = regions.begin(); it != regions.end(); it++)
  {
    it->start = it->start + it->peak - size / 2;
    it->peak = size / 2;
    it->end = it->start + size;
  }
}

void NarrowPeak::mergeOverlappingRegions(bool each_dataset_separately)
{
  if (each_dataset_separately)
  {
    // join overlapping regions within the same dataset
    for (list<Region>::iterator it = regions.begin(); it != regions.end(); it++)
    {
      list<Region>::iterator st = it;
      st++;
      while (st != regions.end() && st->chrom_id == it->chrom_id && st->start <= it->end)
      {
        if (it->dataset_id == st->dataset_id)
        {
          it->end = max(it->end, st->end);
          regions.erase(st++);
        }
        else
          st++;
      }
    }
  }
  else
  {
    // join all overlapping regions
    for (list<Region>::iterator it = regions.begin(); it != regions.end(); it++)
    {
      list<Region>::iterator st = it;
      st++;
      while (st != regions.end() && st->chrom_id == it->chrom_id && st->start <= it->end)
      {
        it->end = max(it->end, st->end);
        regions.erase(st++);
      }
    }
  }
}
