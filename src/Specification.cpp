/*
    Specification.cpp

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
#include <wordexp.h>

#include "Specification.h"

const char *spec_prefix = ".spec";

void Specification::insert_words(char *value, set<string> &set)
{
  char word[1024];
  int n;

  while (sscanf(value, "%s%n", word, &n) > 0)
  {
    set.insert(word);
    value += n;
  }
}

void Specification::push_back_wordexp(char *value, vector<string> &vec)
{
  wordexp_t p;
  if (!wordexp(value, &p, 0))
  {
    for (unsigned int i = 0; i < p.we_wordc; i++)
      vec.push_back(p.we_wordv[i]);
    wordfree(&p);
  }
  else
    vec.push_back(value);
}

void Specification::scan_int(char *value, int &where)
{
  if (sscanf(value, "%d", &where) < 1)
  {
    cerr << "Specification: malformed integer value in line " << lineNum << endl;
    exit(1);
  }
}

void Specification::scan_double(char *value, double &where)
{
  if (sscanf(value, "%lf", &where) < 1)
  {
    cerr << "Specification: malformed floating point value in line " << lineNum << endl;
    exit(1);
  }
}

void Specification::scan_bool(char *value, bool &where)
{
  if (strcasecmp(value, "True") == 0 || strcasecmp(value, "Yes") == 0 || strcasecmp(value, "1") == 0)
    where = true;
  else if (strcasecmp(value, "False") == 0 || strcasecmp(value, "No") == 0 || strcasecmp(value, "0") == 0)
    where = false;
  else
  {
    cerr << "Specification: malformed boolean value in line " << lineNum << endl;
    exit(1);
  }
}

void Specification::scan_RegionMasking(char *value, int &where)
{
  if (strcasecmp(value, "None") == 0)
    where = REGION_MASKING_NONE;
  else if (strcasecmp(value, "Peak") == 0)
    where = REGION_MASKING_PEAK;
  else if (strcasecmp(value, "Majority") == 0)
    where = REGION_MASKING_MAJORITY;
  else
    cout << "Specification: ignoring unknown value \"" << value << "\" in line " << lineNum << endl;
}

void Specification::scan_DatasetList(char *value, vector<pair<string, string> > &where)
{
  FILE *f = fopen(value, "r");
  if (f == NULL)
  {
    cerr << "Could not open file \"" << value << "\" for input" << endl;
    exit(1);
  }

  char line[1024], dataset[1024], fname[1024];
  wordexp_t p;

  while (fgets(line, 1024, f))
    if (sscanf(line, "%s %[^\n]", dataset, fname) == 2)
    {
      if (!wordexp(fname, &p, 0))
      {
        for (unsigned int i = 0; i < p.we_wordc; i++)
          where.push_back(make_pair(dataset, p.we_wordv[i]));
        wordfree(&p);
      }
      else
        where.push_back(make_pair(dataset, fname));
    }
    else
    {
      if (!wordexp(line, &p, 0))
      {
        for (unsigned int i = 0; i < p.we_wordc; i++)
          where.push_back(make_pair(p.we_wordv[i], p.we_wordv[i]));
        wordfree(&p);
      }
      else
        where.push_back(make_pair(line, line));
    }

  fclose(f);
}

void Specification::scan_Dataset(char *value, vector<pair<string, string> > &where)
{
  char dataset[1024], fname[1024];
  wordexp_t p;

  if (sscanf(value, "%s %[^\n]", dataset, fname) == 2)
    {
      if (!wordexp(fname, &p, 0))
      {
        for (unsigned int i = 0; i < p.we_wordc; i++)
          where.push_back(make_pair(dataset, p.we_wordv[i]));
        wordfree(&p);
      }
      else
        where.push_back(make_pair(dataset, fname));
    }
  else
    {
      if (!wordexp(value, &p, 0))
      {
        for (unsigned int i = 0; i < p.we_wordc; i++)
          where.push_back(make_pair(p.we_wordv[i], p.we_wordv[i]));
        wordfree(&p);
      }
      else
        where.push_back(make_pair(value, value));
    }
}

void Specification::scan_output_range(char *value, int &where)
{
  if (strcasecmp(value, "None") == 0)
    where = OUTPUT_RANGE_NONE;
  else if (strcasecmp(value, "Signature") == 0)
    where = OUTPUT_RANGE_SIGNATURE;
  else if (strcasecmp(value, "All") == 0)
    where = OUTPUT_RANGE_ALL;
  else
    cout << "Specification: ignoring unknown value \"" << value << "\" in line " << lineNum << endl;
}

Specification::Specification(char *specFile)
{
  vector<string> context;
  cout << "Reading specification file: " << specFile << endl;

  FILE *f = fopen(specFile, "r");
  if (f == NULL)
  {
    cerr << specFile << ": " << strerror(errno) << endl;
    exit(1);
  }

  // default values of various options
  Options.NumberOfThreads = 1;

  Options.MinMotifInformationContribution = 6.0;
  Options.MaxOverlappingInformationContent = 2.0;
  Options.MaxMotifSpacing = 50;
  Options.ConsiderOrientationsSeparately = true;
  Options.ConsiderMostSignificantComplexOnly = false;

  Options.TargetInstancesThreshold = 100;
  Options.FoldChangeThreshold = 1;
  Options.PValueThreshold = 0.05;

  Options.DimerMotifFlanks = 5;
  Options.ClusteringAcrossDatasets = true;
  Options.ClusteringDistanceConstant = 0;
  Options.ClusteringDistanceMultiplier = 0.15;
  Options.ClusteringOverlapThreshold = 0.2;

  Options.OutputDetailedStats = OUTPUT_RANGE_ALL;
  Options.OutputDimerMotifs = OUTPUT_RANGE_ALL;
  Options.OutputGenomicLocations = OUTPUT_RANGE_ALL;
  Options.GenomicLocationsMaxSpacingDeviation = 0;
  Options.OutputPValueDistribution = true;

  // default output prefix is the specfile prefix (after truncating ".spec") or jutst specfile name
  OutputPrefix = specFile;
  if (strcasecmp(specFile + strlen(specFile) - strlen(spec_prefix), spec_prefix) == 0)
    OutputPrefix.erase(OutputPrefix.size() - strlen(spec_prefix), strlen(spec_prefix));

  char line[1024], key[1024], value[1024];
  lineNum = 0;

  while (fgets(line, 1024, f))
  {
    lineNum++;

    // truncate comments
    char *p = strchr(line, '#');
    if (p) *p = '\0';

    // check for opening/closing tag
    if (sscanf(line, " <%[^>]>", key) == 1)
    {
      if (key[0] == '/')
      {
        if (!context.empty() && string(&key[1]) == context.back())
        {
          context.pop_back();
          continue;
        }
        else
        {
          cerr << "Specification: unexpected closing tag <" << key << "> in line " << lineNum << endl;
          exit(1);
        }
      }
      else
      {
        context.push_back(key);

        if (context.size() == 1 && context.front() == "Genome")
          continue;

        // default values for dataset normalization
        if (context.size() == 1 && context.front() == "StronglySpecificDatasets")
        {
          StronglySpecificDatasets.push_back(DatasetsSection());
          StronglySpecificDatasets.back().RegionSize = 0;
          StronglySpecificDatasets.back().RegionMasking = REGION_MASKING_NONE;
          StronglySpecificDatasets.back().RegionCount = 0;
          continue;
        }

        if (context.size() == 1 && context.front() == "WeaklySpecificDatasets")
        {
          WeaklySpecificDatasets.push_back(DatasetsSection());
          WeaklySpecificDatasets.back().RegionSize = 0;
          WeaklySpecificDatasets.back().RegionMasking = REGION_MASKING_NONE;
          WeaklySpecificDatasets.back().RegionCount = 0;
          continue;
        }

        if (context.size() == 1 && context.front() == "Motifs")
        {
          Motifs.push_back(MotifsSection());
          Motifs.back().Sensitivity = 0.9;
          continue;
        }

        if (context.size() == 1 && context.front() == "Scope")
        {
          Scope.push_back(ScopeSection());
          continue;
        }

        if (context.size() == 1 && context.front() == "Options")
          continue;
      }

      cout << "Specification: ignoring unknown tag <" << context.back() << "> in line " << lineNum << endl;
      continue;
    }

    int count = sscanf(line, " %[^=] = %[^\n]", key, value); // 0 -- empty line, 2 -- key=value pair

    if (count == EOF) continue; // empty line

    if (count != 2)
    {
      cerr << "Specification: incorrect key=value pair in line " << lineNum << endl;
      exit(1);
    }

    // trim whitespace from key
    p = key + strlen(key) - 1;
    while(p > key && isspace(*p)) p--;
    *(p + 1) = '\0';

    if (context.size() == 1 && context.front() == "Genome")
    {
      if (strcasecmp(key, "FastaFile") == 0) { push_back_wordexp(value, Genome.FastaFile); continue; }
      if (strcasecmp(key, "MaskedRegions") == 0) { push_back_wordexp(value, Genome.MaskedRegions); continue; }
    }

    if (context.size() == 1 && context.front() == "StronglySpecificDatasets")
    {
      if (strcasecmp(key, "DatasetList") == 0) { scan_DatasetList(value, StronglySpecificDatasets.back().Dataset); continue; }
      if (strcasecmp(key, "Dataset") == 0) { scan_Dataset(value, StronglySpecificDatasets.back().Dataset); continue; }
      if (strcasecmp(key, "RegionSize") == 0) { scan_int(value, StronglySpecificDatasets.back().RegionSize); continue; }
      if (strcasecmp(key, "RegionMasking") == 0) { scan_RegionMasking(value, StronglySpecificDatasets.back().RegionMasking); continue; }
      if (strcasecmp(key, "RegionCount") == 0) { scan_int(value, StronglySpecificDatasets.back().RegionCount); continue; }
    }

    if (context.size() == 1 && context.front() == "WeaklySpecificDatasets")
    {
      if (strcasecmp(key, "DatasetList") == 0) { scan_DatasetList(value, WeaklySpecificDatasets.back().Dataset); continue; }
      if (strcasecmp(key, "Dataset") == 0) { scan_Dataset(value, WeaklySpecificDatasets.back().Dataset); continue; }
      if (strcasecmp(key, "RegionSize") == 0) { scan_int(value, WeaklySpecificDatasets.back().RegionSize); continue; }
      if (strcasecmp(key, "RegionMasking") == 0) { scan_RegionMasking(value, WeaklySpecificDatasets.back().RegionMasking); continue; }
      if (strcasecmp(key, "RegionCount") == 0) { scan_int(value, WeaklySpecificDatasets.back().RegionCount); continue; }
    }

    if (context.size() == 1 && context.front() == "Motifs")
    {
      if (strcasecmp(key, "Database") == 0) { push_back_wordexp(value, Motifs.back().MotifDatabase); continue; }
      if (strcasecmp(key, "DatabaseSubset") == 0) { push_back_wordexp(value, Motifs.back().MotifSubset); continue; }
      if (strcasecmp(key, "Motif") == 0) { scan_Dataset(value, Motifs.back().Motif); continue; }
      if (strcasecmp(key, "Sensitivity") == 0) { scan_double(value, Motifs.back().Sensitivity); continue; }
    }

    if (context.size() == 1 && context.front() == "Scope")
    {
      if (strcasecmp(key, "Motif1") == 0) { insert_words(value, Scope.back().Motif1); continue; }
      if (strcasecmp(key, "Motif2") == 0) { insert_words(value, Scope.back().Motif2); continue; }
      if (strcasecmp(key, "Dataset") == 0) { insert_words(value, Scope.back().Dataset); continue; }
    }

    if (context.size() == 1 && context.front() == "Options")
    {
      if (strcasecmp(key, "OutputPrefix") == 0) { OutputPrefix = value; continue; }
      if (strcasecmp(key, "NumberOfThreads") == 0) { scan_int(value, Options.NumberOfThreads); continue; }

      if (strcasecmp(key, "MinMotifInformationContribution") == 0) { scan_double(value, Options.MinMotifInformationContribution); continue; }
      if (strcasecmp(key, "MaxOverlappingInformationContent") == 0) { scan_double(value, Options.MaxOverlappingInformationContent); continue; }
      if (strcasecmp(key, "MaxMotifSpacing") == 0) { scan_int(value, Options.MaxMotifSpacing); continue; }
      if (strcasecmp(key, "ConsiderOrientationsSeparately") == 0) { scan_bool(value, Options.ConsiderOrientationsSeparately); continue; }
      if (strcasecmp(key, "ConsiderMostSignificantComplexOnly") == 0) { scan_bool(value, Options.ConsiderMostSignificantComplexOnly); continue; }

      if (strcasecmp(key, "TargetInstancesThreshold") == 0) { scan_int(value, Options.TargetInstancesThreshold); continue; }
      if (strcasecmp(key, "FoldChangeThreshold") == 0) { scan_double(value, Options.FoldChangeThreshold); continue; }
      if (strcasecmp(key, "PValueThreshold") == 0) { scan_double(value, Options.PValueThreshold); continue; }

      if (strcasecmp(key, "DimerMotifFlanks") == 0) { scan_int(value, Options.DimerMotifFlanks); continue; }
      if (strcasecmp(key, "ClusteringAcrossDatasets") == 0) { scan_bool(value, Options.ClusteringAcrossDatasets); continue; }
      if (strcasecmp(key, "ClusteringDistanceConstant") == 0) { scan_double(value, Options.ClusteringDistanceConstant); continue; }
      if (strcasecmp(key, "ClusteringDistanceMultiplier") == 0) { scan_double(value, Options.ClusteringDistanceMultiplier); continue; }
      if (strcasecmp(key, "ClusteringOverlapThreshold") == 0) { scan_double(value, Options.ClusteringOverlapThreshold); continue; }

      if (strcasecmp(key, "OutputDetailedStats") == 0) { scan_output_range(value, Options.OutputDetailedStats); continue; }
      if (strcasecmp(key, "OutputDimerMotifs") == 0) { scan_output_range(value, Options.OutputDimerMotifs); continue; }
      if (strcasecmp(key, "OutputGenomicLocations") == 0) { scan_output_range(value, Options.OutputGenomicLocations); continue; }
      if (strcasecmp(key, "GenomicLocationsMaxSpacingDeviation") == 0) { scan_int(value, Options.GenomicLocationsMaxSpacingDeviation); continue; }
      if (strcasecmp(key, "OutputPValueDistribution") == 0) { scan_bool(value, Options.OutputPValueDistribution); continue; }
    }

    cout << "Specification: ignoring unknown key \"" << key << "\" in line " << lineNum << endl;
  }

  if (!context.empty())
  {
    cerr << "Specification: end of file encountered, but tag <" << context.back() << "> not closed" << endl;
    exit(1);
  }

  // if no Scope is specified, provide the default, empty one (i.e. test all possible combinations)
  if (Scope.empty()) Scope.push_back(ScopeSection());

  cout << "Output will be written to: " << OutputPrefix << ".*" << endl;
  fclose(f);
}
