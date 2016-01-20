/*
    TACO -- Transcription factor Association from Complex Overrepresentation
    (formely: Transcription factor Association from Chromatin Openness)

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
#include <iostream>
#include <math.h>
#include <queue>
#include <sstream>
#include <stdio.h>

#include "Genome.h"
#include "HypothesesSet.h"
#include "NarrowPeak.h"
#include "PositionWeightMatrix.h"
#include "Specification.h"
#include "StatsSet.h"

using namespace std;

const string lead = "TACO -- Transcription factor Association from Complex Overrepresentation";
const string version = "Version 1.0 (git-devel)"; // 1.0 was dated 2013-10-08

Genome fa;
double GC_content;
vector<PositionWeightMatrix> motifs;
map<int, NarrowPeak> target_datasets;
NarrowPeak control_dataset;

vector<PositionWeightMatrix>::iterator motifs_it;
queue<pair<PositionWeightMatrix *, PositionWeightMatrix *> > motif_pair_scope;
map<int, const NarrowPeak *> target_datasets_scope;
list<HypothesesSet::Hypothesis>::iterator hypotheses_it;
pthread_mutex_t data_mutex, cout_mutex;

void readDatasets(Specification &spec)
{
  map<int, int> replicate_to_dataset;

  if (!spec.StronglySpecificDatasets.empty())
    cout << endl << "Reading strongly specific input datasets..." << endl;

  for (vector<Specification::DatasetsSection>::iterator it = spec.StronglySpecificDatasets.begin(); it != spec.StronglySpecificDatasets.end(); it++)
  {
    NarrowPeak this_control_dataset;

    for (vector<pair<string, string> >::iterator jt = it->Dataset.begin(); jt != it->Dataset.end(); jt++)
    {
      int replicate_id = Genome::get_anonymous_dataset_id();
      replicate_to_dataset[replicate_id] = Genome::get_dataset_id(jt->first);
      this_control_dataset.readNarrowPeakFile(replicate_id, jt->second);
    }

    if (it->RegionSize > 0)
      this_control_dataset.setRegionSize(it->RegionSize);
    switch (it->RegionMasking)
    {
      case REGION_MASKING_PEAK: this_control_dataset.removeRegionsWithMaskedPeaks(fa); break;
      case REGION_MASKING_MAJORITY: this_control_dataset.removeMostlyMaskedRegions(fa); break;
    }
    if (it->RegionCount > 0)
      this_control_dataset.restrictToTopSignalRegions(it->RegionCount, true); // for each replicate separately

    // merge strongly specific input dataset replicates into datasets
    for (list<NarrowPeak::Region>::iterator it = this_control_dataset.regions.begin(); it != this_control_dataset.regions.end(); it++)
      it->dataset_id = replicate_to_dataset[it->dataset_id];

    // add the strongly specific input dataset to the control dataset
    for (list<NarrowPeak::Region>::iterator jt = this_control_dataset.regions.begin(); jt != this_control_dataset.regions.end(); jt++)
      control_dataset.regions.push_back(*jt);
  }

  control_dataset.regions.sort();
  control_dataset.mergeOverlappingRegions(true); // for each dataset separately

  if (!spec.StronglySpecificDatasets.empty())
    cout << endl << "Identifying strongly specific target datasets..." << endl;

  for (list<NarrowPeak::Region>::iterator it = control_dataset.regions.begin(); it != control_dataset.regions.end(); it++)
  {
    int new_end = it->end;
    list<NarrowPeak::Region>::iterator st = it;
    st++;

    while (st != control_dataset.regions.end() && st->chrom_id == it->chrom_id && st->start < it->end)
    {
      assert(st->dataset_id != it->dataset_id);
      new_end = min(new_end, st->start);

      st->peak = st->start + st->peak - it->end;
      st->start = it->end;
      if (st->start >= st->end) control_dataset.regions.erase(st++); else st++;
    }

    if (it->start < new_end)
    {
      target_datasets[it->dataset_id].regions.push_back(*it);
      target_datasets[it->dataset_id].regions.back().end = new_end;
    }
  }

  int first_weakly_specific_dataset_id = Genome::get_anonymous_dataset_id();

  if (!spec.WeaklySpecificDatasets.empty())
    cout << endl << "Reading weakly specific input datasets..." << endl;

  for (vector<Specification::DatasetsSection>::iterator it = spec.WeaklySpecificDatasets.begin(); it != spec.WeaklySpecificDatasets.end(); it++)
  {
    int first_dataset_id_in_loop = Genome::get_anonymous_dataset_id();

    for (vector<pair<string, string> >::iterator jt = it->Dataset.begin(); jt != it->Dataset.end(); jt++)
    {
      int replicate_id = Genome::get_anonymous_dataset_id();
      int dataset_id = Genome::get_dataset_id(jt->first);
      if (dataset_id < first_weakly_specific_dataset_id)
      {
        cerr << "Dataset \"" << jt->first << "\" cannot be both in StronglySpecificDatasets and WeaklySpecificDatasets" << endl;
        exit(1);
      }

      replicate_to_dataset[replicate_id] = dataset_id;
      target_datasets[dataset_id].readNarrowPeakFile(replicate_id, jt->second);
    }

    for (map<int, NarrowPeak>::iterator jt = target_datasets.begin(); jt != target_datasets.end(); jt++)
    {
      if (jt->first < first_dataset_id_in_loop) continue; // process only the weakly specific input datasets
      if (it->RegionSize > 0)
        jt->second.setRegionSize(it->RegionSize);
      switch (it->RegionMasking)
      {
        case REGION_MASKING_PEAK: jt->second.removeRegionsWithMaskedPeaks(fa); break;
        case REGION_MASKING_MAJORITY: jt->second.removeMostlyMaskedRegions(fa); break;
      }
      if (it->RegionCount > 0)
        jt->second.restrictToTopSignalRegions(it->RegionCount, true); // for each replicate separately

      // merge weakly specific input dataset replicates into datasets
      for (list<NarrowPeak::Region>::iterator kt = jt->second.regions.begin(); kt != jt->second.regions.end(); kt++)
        kt->dataset_id = replicate_to_dataset[kt->dataset_id];

      jt->second.regions.sort();
      jt->second.mergeOverlappingRegions(true); // for each dataset separately

      // add the weakly specific input dataset to the control dataset
      for (list<NarrowPeak::Region>::iterator kt = jt->second.regions.begin(); kt != jt->second.regions.end(); kt++)
        control_dataset.regions.push_back(*kt);
    }
  }

  // update location_dist statistics for all target datasets
  cout << endl << "Calculating statistics for all target datasets..." << endl;

  for (map<int, NarrowPeak>::iterator it = target_datasets.begin(); it != target_datasets.end(); it++)
  {
    it->second.update_location_dist(fa);
    cout << Genome::dataset_names[it->first] << ": " << it->second.regions.size() << " regions, " << it->second.location_dist[0] << " bp" << endl;

    if (spec.Options.OutputFastaDatasets)
    {
      char fname[1024];
      snprintf(fname, 1024, "%s_%s.fa", spec.OutputPrefix.c_str(), Genome::dataset_names[it->first].c_str());
      it->second.writeFastaFile(fa, fname);
    }
  }

  // sort the control dataset and update its location_dist statistics
  cout << endl << "Calculating statistics for control dataset..." << endl;

  control_dataset.regions.sort();
  control_dataset.mergeOverlappingRegions(false); // merge regions from all target datasets
  control_dataset.update_location_dist(fa);
  cout << "control: " << control_dataset.regions.size() << " regions, " << control_dataset.location_dist[0] << " bp" << endl;

  if (spec.Options.OutputFastaDatasets)
  {
    char fname[1024];
    snprintf(fname, 1024, "%s_control.fa", spec.OutputPrefix.c_str());
    control_dataset.writeFastaFile(fa, fname);
  }

  if (control_dataset.location_dist[0] == 0)
  {
    cerr << "Control dataset cannot be empty" << endl;
    exit(1);
  }

  // calculate the GC-content within control dataset
  GC_content = fa.calculateGC_content(&control_dataset);
  cout << "GC-content: " << GC_content << endl;
}

void readMotifSubset(set<string> &subset, const char *fname)
{
  FILE *f = fopen(fname, "r");
  if (f == NULL)
  {
    cerr << "Could not open file \"" << fname << "\" for input" << endl;
    exit(1);
  }

  char line[1024], word[1024];
  while (fgets(line, 1024, f))
  {
    char *value = line;
    int n;

    while (sscanf(value, "%s%n", word, &n) > 0)
    {
      subset.insert(word);
      value += n;
    }
  }

  fclose(f);
}

void fixMotifAccession(map<string, int> &accessions, PositionWeightMatrix &motif)
{
  if (accessions[motif.accession] > 0)
    motif.accession += '_' + static_cast<ostringstream*>(&(ostringstream() << accessions[motif.accession]))->str();
  accessions[motif.accession]++;

  cout << "Motif " << motif.accession << ": id \"" << motif.id << "\", name \"" << motif.name << "\", " << motif.length << " bp, information content " << motif.inf_content << " bits" << endl;
}

void readMotifDatabase(map<string, int> &accessions, vector<PositionWeightMatrix> &motifs, Specification::MotifsSection &section, set<string> &subset, const char *fname)
{
  FILE *f = fopen(fname, "r");
  if (f == NULL)
  {
    cerr << "Could not open file \"" << fname << "\" for input" << endl;
    exit(1);
  }
  int lineNum = 0;

  try
  {
    while (true)
    {
      PositionWeightMatrix motif = PositionWeightMatrix(GC_content, f, fname, &lineNum, section.Pseudocounts);

      if (!section.MotifSubset.empty() && subset.find(motif.accession) == subset.end()) continue;

      fixMotifAccession(accessions, motif);
      motif.estimateThreshold(section.Sensitivity);
      motifs.push_back(motif);
    }
  }
  catch (PositionWeightMatrix::EndOfFile &e) {};

  fclose(f);
}

void readMotif(map<string, int> &accessions, vector<PositionWeightMatrix> &motifs, Specification::MotifsSection &section, const string accession, const char *fname)
{
  PositionWeightMatrix motif = PositionWeightMatrix(GC_content, fname, accession, section.Pseudocounts);
  fixMotifAccession(accessions, motif);
  motif.estimateThreshold(section.Sensitivity);
  motifs.push_back(motif);
}

void readMotifs(Specification &spec)
{
  map<string, int> accessions;
  cout << endl << "Reading transcription factor motifs..." << endl;

  for (vector<Specification::MotifsSection>::iterator it = spec.Motifs.begin(); it != spec.Motifs.end(); it++)
  {
    set<string> subset;
    for (vector<string>::iterator jt = it->MotifSubset.begin(); jt != it->MotifSubset.end(); jt++)
      readMotifSubset(subset, jt->c_str());

    if (!it->MotifSubset.empty())
      cout << "Motif subset: " << subset.size() << " identifiers" << endl;

    for (vector<string>::iterator jt = it->MotifDatabase.begin(); jt != it->MotifDatabase.end(); jt++)
      readMotifDatabase(accessions, motifs, *it, subset, jt->c_str());

    for (vector<pair<string, string> >::iterator jt = it->Motif.begin(); jt != it->Motif.end(); jt++)
      readMotif(accessions, motifs, *it, jt->first, jt->second.c_str());
  }

  cout << "Total number of motifs: " << motifs.size() << endl;
}

void *motifScanThread(void *arg)
{
  while (true)
  {
    pthread_mutex_lock(&data_mutex);
    if (motifs_it == motifs.end()) { pthread_mutex_unlock(&data_mutex); return NULL; }
    PositionWeightMatrix &par = *motifs_it++;
    pthread_mutex_unlock(&data_mutex);

    par.scan(fa, &control_dataset);
    pthread_mutex_lock(&cout_mutex);
    cout << "Motif " << par.accession << ": " << par.matches.size() << " matches" << endl;
    pthread_mutex_unlock(&cout_mutex);
  }
}

void *experimentThread(void *arg)
{
  HypothesesSet *hs = (HypothesesSet *) arg;

  while (true)
  {
    pthread_mutex_lock(&data_mutex);
    if (motif_pair_scope.empty()) { pthread_mutex_unlock(&data_mutex); return NULL; }
    pair<PositionWeightMatrix *, PositionWeightMatrix *> par = motif_pair_scope.front();
    motif_pair_scope.pop();
    pthread_mutex_unlock(&data_mutex);

    hs->processPair(par.first, par.second, target_datasets_scope, &control_dataset);
  }
}

void *postprocessThreadDimerMotifSimilarity(void *arg)
{
  HypothesesSet *hs = (HypothesesSet *) arg;

  while (true)
  {
    pthread_mutex_lock(&data_mutex);
    if (hypotheses_it == hs->hypotheses_end()) { pthread_mutex_unlock(&data_mutex); return NULL; }
    HypothesesSet::Hypothesis &par = *hypotheses_it++;
    pthread_mutex_unlock(&data_mutex);

    par.calculateDimerMotifSimilarity(fa, motifs, GC_content, hs->spec.Options.DimerMotifFlanks);
  }
}

void *postprocessThreadDimerMotif(void *arg)
{
  HypothesesSet *hs = (HypothesesSet *) arg;

  while (true)
  {
    pthread_mutex_lock(&data_mutex);
    if (hypotheses_it == hs->hypotheses_end()) { pthread_mutex_unlock(&data_mutex); return NULL; }
    HypothesesSet::Hypothesis &par = *hypotheses_it++;
    pthread_mutex_unlock(&data_mutex);

    par.calculateDimerMotif(fa, GC_content, hs->spec.Options.DimerMotifFlanks);
  }
}

int main(int argc, char *argv[])
{
  if (argc < 2 || argv[1][0] == '-')
  {
    cerr << lead << endl << version << endl << endl << "See http://bioputer.mimuw.edu.pl/taco/ for license and documentation." << endl << endl << "Usage:  " << argv[0] << " <specification_file>" << endl << endl;
    exit(1);
  }

  cout << lead << endl << version << endl << endl;
  Specification spec(argv[1]);
  cout << endl << "Reading genome sequence..." << endl;
  HypothesesSet hs(spec);
  char fname[1024];

  for (vector<string>::iterator it = spec.Genome.FastaFile.begin(); it != spec.Genome.FastaFile.end(); it++)
    fa.readFasta(*it);
  for (vector<string>::iterator it = spec.Genome.MaskedRegions.begin(); it != spec.Genome.MaskedRegions.end(); it++)
    fa.applyMask(*it);

  readDatasets(spec);
  readMotifs(spec);

  pthread_mutex_init(&data_mutex, NULL);
  pthread_mutex_init(&cout_mutex, NULL);
  pthread_t *pth = new pthread_t[spec.Options.NumberOfThreads];

  // scan the genome limited to control_dataset for motif occurrences
  cout << endl << "Scanning the control dataset for individual motif occurrences..." << endl;
  motifs_it = motifs.begin();
  for (int i = 0; i < spec.Options.NumberOfThreads; i++)
    pthread_create(&pth[i], NULL, motifScanThread, NULL);
  for (int i = 0; i < spec.Options.NumberOfThreads; i++)
    pthread_join(pth[i], NULL);

  cout << endl << "Predicting cell-type-specific overrepresented motif complexes..." << endl;
  for (vector<Specification::ScopeSection>::iterator it = spec.Scope.begin(); it != spec.Scope.end(); it++)
  {
    // insert all combinations of (M1, M2) into the parameter queue
    for (vector<PositionWeightMatrix>::iterator jt = motifs.begin(); jt != motifs.end(); jt++)
      for (vector<PositionWeightMatrix>::iterator kt = jt; kt != motifs.end(); kt++)
      {
        bool jt_scope1 = it->Motif1.empty() || it->Motif1.find(jt->accession) != it->Motif1.end();
        bool kt_scope2 = it->Motif2.empty() || it->Motif2.find(kt->accession) != it->Motif2.end();
        if (jt_scope1 && kt_scope2)
        {
          motif_pair_scope.push(make_pair(&*jt, &*kt));
          continue;
        }

        bool kt_scope1 = it->Motif1.empty() || it->Motif1.find(kt->accession) != it->Motif1.end();
        bool jt_scope2 = it->Motif2.empty() || it->Motif2.find(jt->accession) != it->Motif2.end();
        if (kt_scope1 && jt_scope2)
          motif_pair_scope.push(make_pair(&*kt, &*jt));
      }

    // set the scope of datasets considered
    target_datasets_scope.clear();
    for (map<int, NarrowPeak>::iterator jt = target_datasets.begin(); jt != target_datasets.end(); jt++)
      if (it->Dataset.empty() || it->Dataset.find(Genome::dataset_names[jt->first]) != it->Dataset.end())
        target_datasets_scope[jt->first] = &jt->second;

    // main experiment processing: consider all the combinations of (M1, M2, dataset)
    for (int i = 0; i < spec.Options.NumberOfThreads; i++)
      pthread_create(&pth[i], NULL, experimentThread, &hs);
    for (int i = 0; i < spec.Options.NumberOfThreads; i++)
      pthread_join(pth[i], NULL);
  }

  cout << "Number of hypotheses considered: " << hs.getNumberOfHypotheses() << endl;
  // remove insignificant hypotheses, then save the predictions and their genomic locations
  hs.removeInsignificantHypotheses();
  cout << "Number of accepted hypotheses: " << hs.getNumberOfAcceptedHypotheses() << endl;

  if (spec.Options.ClusteringDistanceConstant > 0.0 || spec.Options.ClusteringDistanceMultiplier > 0.0 || spec.Options.AnnotateDimerMotifSimilarity)
  {
    // before clustering, calculate the dimer motifs and possibly find similar individual motifs
    hypotheses_it = hs.hypotheses_begin();

    if (spec.Options.AnnotateDimerMotifSimilarity)
    {
      cout << endl << "Calculating dimer motifs and annotating their similarity to individual motifs..." << endl;
      for (int i = 0; i < spec.Options.NumberOfThreads; i++)
        pthread_create(&pth[i], NULL, postprocessThreadDimerMotifSimilarity, &hs);
    }
    else
    {
      cout << endl << "Calculating dimer motifs..." << endl;
      for (int i = 0; i < spec.Options.NumberOfThreads; i++)
        pthread_create(&pth[i], NULL, postprocessThreadDimerMotif, &hs);
    }

    for (int i = 0; i < spec.Options.NumberOfThreads; i++)
      pthread_join(pth[i], NULL);
  }

  cout << endl << "Clustering the overrepresented motif complexes..." << endl;
  // cluster the hypotheses, then save the clustering results and dimer motifs
  hs.clusterHypotheses();
  cout << "Number of distinct predictions: " << hs.getNumberOfClusters() << endl;

  cout << endl << "Writing output files..." << endl;
  snprintf(fname, 1024, "%s.tab", spec.OutputPrefix.c_str());
  hs.writeClusteringResultsFile(fname);
  snprintf(fname, 1024, "%s.stats", spec.OutputPrefix.c_str());
  hs.writeDetailedStatsFile(fname);
  snprintf(fname, 1024, "%s.pwms", spec.OutputPrefix.c_str());
  hs.writeDimerMotifsFile(fname);
  snprintf(fname, 1024, "%s.hits", spec.OutputPrefix.c_str());
  hs.writeGenomicLocationsFile(fname);
  snprintf(fname, 1024, "%s.pval", spec.OutputPrefix.c_str());
  hs.writePValueDistributionFile(fname);

  delete[] pth;
  pthread_mutex_destroy(&cout_mutex);
  pthread_mutex_destroy(&data_mutex);

  cout << endl << "Done! Thank you for using TACO." << endl;
  return 0;
}
