/*
    Specification.h

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

#ifndef _SPECIFICATION_H
#define _SPECIFICATION_H

#include <set>
#include <string>

#include "Genome.h"

using namespace std;

#define REGION_MASKING_NONE 0
#define REGION_MASKING_PEAK 1
#define REGION_MASKING_MAJORITY 2

#define OUTPUT_RANGE_NONE 0
#define OUTPUT_RANGE_SIGNATURE 1
#define OUTPUT_RANGE_ALL 2

class Specification
{
  public:
    Specification(char *specFile);

    int lineNum;
    string OutputPrefix;

    struct
    {
      vector<string> FastaFile, MaskedRegions;
    } Genome;

    struct DatasetsSection
    {
      vector<pair<string, string> > Dataset;
      int RegionSize, RegionMasking, RegionCount;
    };
    vector<DatasetsSection> StronglySpecificDatasets, WeaklySpecificDatasets;

    struct MotifsSection
    {
      vector<string> MotifDatabase, MotifSubset;
      vector<pair<string, string> > Motif;
      double Sensitivity;
    };
    vector<MotifsSection> Motifs;

    struct ScopeSection
    {
      set<string> Motif1, Motif2, Dataset;
    };
    vector<ScopeSection> Scope;

    struct
    {
      int NumberOfThreads;

      double MinMotifInformationContribution, MaxOverlappingInformationContent;
      int MaxMotifSpacing;
      bool ConsiderOrientationsSeparately, ConsiderMostSignificantComplexOnly;

      int TargetInstancesThreshold;
      double FoldChangeThreshold, PValueThreshold;

      int DimerMotifFlanks;
      bool ClusteringAcrossDatasets;
      double ClusteringDistanceConstant, ClusteringDistanceMultiplier, ClusteringOverlapThreshold;

      int OutputDetailedStats, OutputDimerMotifs, OutputGenomicLocations;
      int GenomicLocationsMaxSpacingDeviation;
      bool OutputPValueDistribution;
    } Options;

  private:
    void insert_words(char *value, set<string> &map);
    void push_back_wordexp(char *value, vector<string> &vec);
    void scan_int(char *value, int &where);
    void scan_double(char *value, double &where);
    void scan_bool(char *value, bool &where);
    void scan_RegionMasking(char *value, int &where);
    void scan_DatasetList(char *value, vector<pair<string, string> > &where);
    void scan_Dataset(char *value, vector<pair<string, string> > &where);
    void scan_output_range(char *value, int &where);
};

#endif
