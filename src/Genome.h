/*
    Genome.h

    Copyright (C) 2011-2014  Aleksander Jankowski <ajank@mimuw.edu.pl>

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

#ifndef _GENOME_H
#define _GENOME_H

#include <map>
#include <string>
#include <vector>

using namespace std;

class NarrowPeak;

class Genome
{
  public:
    static int get_dataset_id(const string &dataset);
    static int get_anonymous_dataset_id();
    static int get_chrom_id(const string &chrom);
    void readFasta(const string &filename);
    void applyMask(const string &filename);
    inline static bool isUnmasked(char n)
    {
      switch (n)
      {
        case 'A': return true;
        case 'C': return true;
        case 'G': return true;
        case 'T': return true;
        default: return false;
      }
    }
    inline static bool isUnmaskedOrSingletonN(const vector<char> &region, int i)
    {
      return isUnmasked(region[i]) || (i > 0 && isUnmasked(region[i - 1]) && (unsigned int) i + 1 < region.size() && isUnmasked(region[i + 1]));
    }
    inline static int nucleotideToIndex(char n) // with masking
    {
      switch (n)
      {
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        default: return 0;
      }
    }
    inline static char indexToNucleotide(int i)
    {
      switch (i)
      {
        case 1: return 'A';
        case 2: return 'C';
        case 3: return 'G';
        case 4: return 'T';
        default: return 'N';
      }
    }
    inline static int complementIndex(int i)
    {
      if (i == 0) return i; else return 5 - i;
    }
    inline static char complementNucleotide(char n)
    {
      switch (n)
      {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
      }
    }

    map<int, vector<char> > chroms;
    double calculateGC_content(NarrowPeak *scope = NULL) const;
    // preferably should be read-only:
    static map<string, int> chrom_ids;
    static map<int, string> chrom_names;
    static map<string, int> dataset_ids;
    static map<int, string> dataset_names;

  private:
    static int chrom_id_count;
    static int dataset_id_count;
};

#endif
