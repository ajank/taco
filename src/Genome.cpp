/*
    Genome.cpp

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

#include <cstdlib>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include "Genome.h"
#include "NarrowPeak.h"

#define nti(nucleotide) nucleotideToIndex(nucleotide)

map<string, int> Genome::chrom_ids;
map<int, string> Genome::chrom_names;
int Genome::chrom_id_count = 0;

map<string, int> Genome::dataset_ids;
map<int, string> Genome::dataset_names;
int Genome::dataset_id_count = 0;

int Genome::get_dataset_id(const string &dataset)
{
  int dataset_id;
  map<string, int>::iterator it = dataset_ids.find(dataset);

  if (it == dataset_ids.end()) dataset_names[dataset_id = dataset_ids[dataset] = ++dataset_id_count] = dataset;
  else dataset_id = it->second;

  return dataset_id;
}

int Genome::get_anonymous_dataset_id()
{
  return ++dataset_id_count;
}

int Genome::get_chrom_id(const string &chrom)
{
  int chrom_id;
  map<string, int>::iterator it = chrom_ids.find(chrom);

  if (it == chrom_ids.end()) chrom_names[chrom_id = chrom_ids[chrom] = ++chrom_id_count] = chrom;
  else chrom_id = it->second;

  return chrom_id;
}

void Genome::readFasta(const string &filename)
{
  cout << "Reading FASTA file: " << filename << endl;
  FILE *f = fopen(filename.c_str(), "r");
  if (f == NULL)
  {
    cerr << filename << ": " << strerror(errno) << endl;
    exit(1);
  }

  char line[1024], chrom_name[10240];
  vector <char>* this_chrom = NULL;

  while (fgets(line, 1024, f))
  {
    if (sscanf(line, "> %s", chrom_name) == 1)
    {
      if (this_chrom != NULL) cout << this_chrom->size() << " bp" << endl;
      cout << chrom_name << ": ";
      this_chrom = &chroms[get_chrom_id(chrom_name)];
    }
    else if (this_chrom != NULL)
    {
      for (char *nucl = line; *nucl != '\0'; nucl++)
        if (isalpha(*nucl)) this_chrom->push_back(*nucl);
    }
  }
  if (this_chrom != NULL) cout << this_chrom->size() << " bp" << endl;
  fclose(f);
}

void Genome::applyMask(const string &filename)
{
  cout << "Applying region mask: " << filename << endl;
  FILE *f = fopen(filename.c_str(), "r");
  if (f == NULL)
  {
    cerr << filename << ": " << strerror(errno) << endl;
    exit(1);
  }

  char line[10240], chrom_name[1024];
  int chromStart, chromEnd, lineNum = 0, numMasked = 0;

  while (fgets(line, 10240, f))
  {
    lineNum++;
    if (sscanf(line, "%s%d%d", chrom_name, &chromStart, &chromEnd) != 3)
      cerr << filename << ": cannot parse line " << lineNum << endl;
    else
    {
      map<string, int>::iterator it = chrom_ids.find(chrom_name);
      if (it != chrom_ids.end())
      {
        map<int, vector<char> >::iterator jt = chroms.find(it->second);
        if (jt != chroms.end())
        {
          numMasked++;
          chromStart = max(0, chromStart);
          chromEnd = min((int) jt->second.size(), chromEnd);

          for (int k = chromStart; k < chromEnd; k++)
            jt->second[k] = tolower(jt->second[k]);
        }
      }
    }
  }

  cout << "Number of masked regions: " << numMasked << endl;
  fclose(f);
}

double Genome::calculateGC_content(NarrowPeak *scope) const
{
  unsigned int counts[5];
  for (int i = 0; i < 5; i++)
    counts[i] = 0;

  if (scope)
    for (list<NarrowPeak::Region>::const_iterator jt = scope->regions.begin(); jt != scope->regions.end(); jt++)
    {
      map<int, vector<char> >::const_iterator it = chroms.find(jt->chrom_id);
      if (it == chroms.end())
      {
        cerr << "Genome: chromosome \"" << chrom_names[jt->chrom_id] << "\" not found in the genome" << endl;
        exit(1);
      }

      for (int j = jt->start; j < jt->end; j++)
        counts[nti(it->second[j])]++;
    }
  else
    for (map<int, vector<char> >::const_iterator it = chroms.begin(); it != chroms.end(); it++)
      for (vector<char>::const_iterator jt = it->second.begin(); jt != it->second.end(); jt++)
        counts[nti(*jt)]++;

  return (double) (counts[nti('C')] + counts[nti('G')]) / (counts[nti('A')] + counts[nti('C')] + counts[nti('G')] + counts[nti('T')]);
}
