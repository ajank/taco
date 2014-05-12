/*
    HypothesesSet.cpp

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
#include <errno.h>
#include <fstream>
#include <iostream>
#include <Rmath.h>

#include "HypothesesSet.h"
#include "StatsSet.h"

const double log10_raw_p_value_sampling = 0.05;

#define nti(nucleotide) Genome::nucleotideToIndex(nucleotide)
#define ots(same_orientation) (same_orientation ? "forward" : "reverse")

bool HypothesesSet::Hypothesis::operator<(const Hypothesis &other) const
{
  if (this->cluster_id < other.cluster_id) return true;
  else if (this->cluster_id > other.cluster_id) return false;
  else return this->raw_log_p_value < other.raw_log_p_value;
}

bool HypothesesSet::Hypothesis::ptr_compare(const Hypothesis *lhs, const Hypothesis *rhs)
{
  return *lhs < *rhs;
}

void HypothesesSet::Hypothesis::calculateStructure()
{
  // at each of the overlapping positions, choose either M1 or M2 column to withdraw at this position
  overlap_inf_content = 0.0; // minimal sum of IC of withdrawn motif matrix columns
  M1_inf_contribution = 0.0;
  M2_inf_contribution = 0.0;
  pair_start = min(0, offset);
  pair_end = max(M1->length, M2->length + offset);

  double M1_mass_center = 0.0, M2_mass_center = 0.0;

  for (int i = pair_start; i < pair_end; i++)
  {
    double M1_pIC = 0.0, M2_pIC = 0.0;

    if (i >= 0 && i < M1->length)
      M1_pIC = M1->position_inf_content[i];

    if (i >= offset && i < M2->length + offset)
    {
      if (same_orientation)
        M2_pIC = M2->position_inf_content[i - offset];
      else
        M2_pIC = M2->position_inf_content[M2->length - 1 - i + offset];
    }

    if (M1_pIC > M2_pIC)
    {
      M1_mass_center += i * M1_pIC;
      M1_inf_contribution += M1_pIC;
      overlap_inf_content += M2_pIC;
    }
    else if (M1_pIC < M2_pIC)
    {
      overlap_inf_content += M1_pIC;
      M2_mass_center += i * M2_pIC;
      M2_inf_contribution += M2_pIC;
    }
    else
    {
      M1_mass_center += i * M1_pIC / 2;
      M1_inf_contribution += M1_pIC / 2;
      overlap_inf_content += (M1_pIC + M2_pIC) / 2;
      M2_mass_center += i * M2_pIC / 2;
      M2_inf_contribution += M2_pIC / 2;
    }
  }

  M1_mass_center /= M1_inf_contribution;
  M2_mass_center /= M2_inf_contribution;

  M1_goes_first = M1_mass_center <= M2_mass_center;
}

void HypothesesSet::Hypothesis::calculateDimerMotif(const Genome &fa, double GC_content, int dimer_motif_margin)
{
  dimer = PositionWeightMatrix(GC_content);
  vector<string> seqs;

  for (vector<PositionWeightMatrix::MotifMatch>::iterator jt = M1_paired_matches.begin(); jt != M1_paired_matches.end(); jt++)
  {
    map<int, vector<char> >::const_iterator it = fa.chroms.find(jt->chrom_id);
    if (it == fa.chroms.end()) continue;
    const vector<char> &chrom = it->second;

    int start = get_start(*jt);
    int actual_start = max(0, start - dimer_motif_margin);
    int end = get_end(*jt);
    int actual_end = min((int) chrom.size(), end + dimer_motif_margin);
    string seq;

    if (jt->strand == '+')
    {
      for (int i = start - dimer_motif_margin; i < actual_start; i++)
        seq.push_back('N');
      for (int i = actual_start; i < actual_end; i++)
        seq.push_back(toupper(chrom[i])); // toupper leaves all the masking (of coding regions, repetitive sequences etc.) out
      for (int i = actual_end; i < end + dimer_motif_margin; i++)
        seq.push_back('N');
    }
    else
    {
      for (int i = end + dimer_motif_margin - 1; i > actual_end - 1; i--)
        seq.push_back('N');
      for (int i = actual_end - 1; i > actual_start - 1; i--)
        seq.push_back(Genome::complementNucleotide(toupper(chrom[i])));
      for (int i = actual_start - 1; i > start - dimer_motif_margin - 1; i--)
        seq.push_back('N');
    }

    seqs.push_back(seq);
  }

  for (int i = 0; i < get_length() + 2 * dimer_motif_margin; i++)
  {
    PositionWeightMatrix::NucleotideValues fd;
    for (int j = 0; j < 5; j++)
      fd.values[j] = 0;

    for (vector<string>::iterator jt = seqs.begin(); jt != seqs.end(); jt++)
      fd.values[Genome::nti((*jt)[i])]++;

    dimer.addPosition(fd);
  }
}

// before clustering, calculate the emprirical motif and find similar individual motifs
void HypothesesSet::Hypothesis::calculateDimerMotifSimilarity(const Genome &fa, const vector<PositionWeightMatrix> &motifs, double GC_content, int dimer_motif_margin)
{
  if (M1_paired_matches.empty()) return; // no genomic instances, no dimer motif

  // calculate the dimer motif
  calculateDimerMotif(fa, GC_content, dimer_motif_margin);

  // compare the motif pair with all the single motifs
  for (vector<PositionWeightMatrix>::const_iterator jt = motifs.begin(); jt != motifs.end(); jt++)
  {
    bool same_orientation = true;
    while (true) // iterating same_orientation over {true, false}
    {
      for (int offset = -jt->length + dimer_motif_margin + 1; offset < dimer.length - 2 * dimer_motif_margin; offset++)
      {
        double dist = dimer.EuclideanDistanceSquared(*jt, offset, same_orientation);
        if (dist < similarity.distance)
        {
          similarity.motif = &*jt;
          similarity.offset = offset - dimer_motif_margin;
          similarity.same_orientation = same_orientation;
          similarity.distance = dist;
        }
      }
      if (same_orientation) same_orientation = false; else break;
    }
  }
}

HypothesesSet::HypothesesSet(const Specification &spec) : spec(spec)
{
  num_hypotheses = 0;
  hypotheses_size = 0;
  hypotheses_previous_size = 0;
  max_log_p_value = log(spec.Options.PValueThreshold);
  num_clusters = 0;
  pthread_mutex_init(&mutex, NULL);
  pthread_mutex_init(&remove_mutex, NULL);
}

void HypothesesSet::processPair(PositionWeightMatrix *M1, PositionWeightMatrix *M2, const map<int, const NarrowPeak *> &target_datasets, const NarrowPeak *control_dataset)
{
  vector<vector<bool> > orientation_sets;

  if (spec.Options.ConsiderOrientationsSeparately)
  {
    orientation_sets.push_back(vector<bool> (1, true));
    orientation_sets.push_back(vector<bool> (1, false));
  }
  else
  {
    vector<bool> orientations;
    orientations.push_back(true);
    orientations.push_back(false);
    orientation_sets.push_back(orientations);
  }

  Hypothesis h;
  h.M1 = M1;
  h.M2 = M2;
  h.cluster_id = 0;
  h.clustering_status = STATUS_UNRESOLVED;
  h.similarity.distance = INFINITY;
  h.raw_log_p_value = INFINITY;

  // return early if the information content of any of the motifs is less than the minimal information contribution
  if (M1->inf_content < spec.Options.MinMotifInformationContribution || M2->inf_content < spec.Options.MinMotifInformationContribution) return;

  StatsSet ss_control(M1, M2, spec.Options.MaxMotifSpacing);
  ss_control.addMotifMatchesFromRegions(control_dataset);
  ss_control.calculateStats();

  for (map<int, const NarrowPeak *>::const_iterator it = target_datasets.begin(); it != target_datasets.end(); it++)
  {
    h.dataset_id = it->first;
    const NarrowPeak *target_dataset = it->second;
    StatsSet ss_target(M1, M2, spec.Options.MaxMotifSpacing);
    ss_target.addMotifMatchesFromRegions(target_dataset);
    ss_target.calculateStats();

    list<Hypothesis> hl; // gathers all the hypotheses for (M1, M2, dataset)
    double min_raw_log_p_value = INFINITY;
    list<Hypothesis>::iterator min_raw_log_p_value_hypothesis = hl.end();

    for (vector<vector<bool> >::iterator orientations = orientation_sets.begin(); orientations != orientation_sets.end(); orientations++)
    {
      list<Hypothesis> ohl; // gathers all the hypotheses for (M1, M2, dataset, orientations)
      // if ConsiderOrientationsSeparately = True, orientations can be {same} or {opposite}
      // otherwise, orientations are {same, opposite} and hl will eventually be the same as ohl
      long int sum_target_instances = 0, sum_target_N = 0, sum_control_instances = 0, sum_control_N = 0;

      for (vector<bool>::const_iterator jt = orientations->begin(); jt != orientations->end(); jt++)
      {
        h.same_orientation = *jt;
        for (h.offset = ss_target.min_offset; h.offset <= ss_target.max_offset; h.offset++)
        {
          int size = abs(h.offset + M2->length / 2 - M1->length / 2); // for calculating motif midpoints, cf. PositionWeightMatrix::scan

          Stats *target = ss_target.getStats(h.offset, h.same_orientation);
          Stats *control = ss_control.getStats(h.offset, h.same_orientation);

          h.target_instances = target->getHits();
          h.target_N = target_dataset->location_dist[size];
          h.control_instances = control->getHits();
          h.control_N = control_dataset->location_dist[size];
          h.calculateStructure();

          if (h.overlap_inf_content > spec.Options.MaxOverlappingInformationContent || h.M1_inf_contribution < spec.Options.MinMotifInformationContribution || h.M2_inf_contribution < spec.Options.MinMotifInformationContribution)
            h.hypothesis_id = HYPOTHESIS_NOT_CONSIDERED; // disallowed motif complex structures
          else
          {
            if ((M1 == M2) && h.same_orientation & (h.offset < 0))
              h.hypothesis_id = HYPOTHESIS_NOT_CONSIDERED; // don't consider the redundant hypotheses in homodimer case
            else
              h.hypothesis_id = HYPOTHESIS_CONSIDERED;

            sum_target_instances += h.target_instances;
            sum_target_N += h.target_N;
            sum_control_instances += h.control_instances;
            sum_control_N += h.control_N;
          }

          ohl.push_back(h);
        }
      }

      // now we have all the sums to calculate p-values
      double prob_base = ((double) sum_target_instances / sum_target_N) / ((double) sum_control_instances / sum_control_N);

      for (list<Hypothesis>::iterator jt = ohl.begin(); jt != ohl.end(); jt++)
      {
        jt->prob = ((double) jt->control_instances / jt->control_N) * prob_base;
        jt->fold_change = jt->target_instances / (jt->prob * jt->target_N);

        // don't consider the hypothesis if f_{12} / b_{12} is smaller than a given threshold
        // (default: FrequencyRatioThreshold = 1.0, i.e. exclude motif pairs which are less frequent in the foreground than in the background)
        // or we have encountered a singularity (sum_target_N == 0, sum_control_instances == 0 or control_N == 0)
        if (prob_base < spec.options.FrequencyRatioThreshold || !isfinite(jt->prob)) jt->hypothesis_id = HYPOTHESIS_NOT_CONSIDERED;

        if (jt->hypothesis_id == HYPOTHESIS_CONSIDERED)
        {
          // consider the hypothesis, increasing the hypotheses count
          pthread_mutex_lock(&mutex);
          jt->hypothesis_id = ++num_hypotheses;
          pthread_mutex_unlock(&mutex);

          // calculate raw p-value, i.e. P(target_instances >= X) where X ~ B(target_N, prob)
          jt->raw_log_p_value = pbinom(jt->target_instances - 1, jt->target_N, jt->prob, 0 /* lower.tail = F */, 1 /* log = T */);

          // save the raw p-value distribution for qq-plot
          int addr = (int) (-get_log10_raw_p_value(*jt) / log10_raw_p_value_sampling);
          pthread_mutex_lock(&mutex);
          if ((int) log10_raw_p_value_dist.size() <= addr)
            log10_raw_p_value_dist.resize(addr + 1, 0);
          log10_raw_p_value_dist[addr]++;
          pthread_mutex_unlock(&mutex);

          // reject the hypothesis it if the motif complex is not frequent enough
          if (jt->target_instances < spec.Options.TargetInstancesThreshold || jt->fold_change < spec.Options.FoldChangeThreshold)
            jt->clustering_status = STATUS_REJECTED;
          else if (jt->raw_log_p_value < min_raw_log_p_value)
          {
            min_raw_log_p_value = jt->raw_log_p_value;
            min_raw_log_p_value_hypothesis = jt;
          }
        }
        else // mark not considered hypotheses as rejected
          jt->clustering_status = STATUS_REJECTED;
      }

      hl.splice(hl.end(), ohl);
    }

    // if all the raw log p-values are equal to INFINITY, choose the removal_hypothesis arbitrarily
    if (min_raw_log_p_value_hypothesis == hl.end())
      min_raw_log_p_value_hypothesis = hl.begin();

    pthread_mutex_lock(&mutex);
    if (is_p_value_significant(min_raw_log_p_value))
    {
      pthread_mutex_unlock(&mutex);
      int size = 0;

      for (list<Hypothesis>::iterator jt = hl.begin(); jt != hl.end(); jt++)
      {
        size++;
        // group together the hypotheses for given (M1, M2, dataset)
        jt->removal_raw_log_p_value = min_raw_log_p_value;
        jt->removal_hypothesis = &*min_raw_log_p_value_hypothesis;

        if (spec.Options.ConsiderMostSignificantComplexOnly && jt != min_raw_log_p_value_hypothesis)
          jt->clustering_status = STATUS_REJECTED;

        if (!is_p_value_significant(jt->raw_log_p_value))
          jt->clustering_status = STATUS_REJECTED;

        // and for the ones that look significant for now, save their genomic instances
        if (jt->clustering_status != STATUS_REJECTED)
        {
          ss_target.returnPairedMatches(jt->M1_paired_matches, jt->offset, jt->same_orientation);
          assert((long int) jt->M1_paired_matches.size() == jt->target_instances);

          int spacing_sign = jt->M1_goes_first ? 1 : -1;
          for (int spacing_deviation = 1; spacing_deviation <= spec.Options.GenomicLocationsMaxSpacingDeviation; spacing_deviation++)
            ss_target.returnPairedMatches(jt->M1_spaced_matches[spacing_deviation], jt->offset + spacing_sign * spacing_deviation, jt->same_orientation);
        }
      }

      pthread_mutex_lock(&mutex);
      hypotheses.splice(hypotheses.end(), hl);
      hypotheses_size += size;
    }
    pthread_mutex_unlock(&mutex);
  }

  pthread_mutex_lock(&mutex);
  if (hypotheses_size > 2 * hypotheses_previous_size)
  {
    pthread_mutex_unlock(&mutex);
    removeInsignificantHypotheses();
  }
  else
    pthread_mutex_unlock(&mutex);
}

void HypothesesSet::removeInsignificantHypotheses()
{
  // recalculate p-values and remove insignificant hypotheses
  if (pthread_mutex_trylock(&remove_mutex) == EBUSY) return; // resign if another thread is doing it
  pthread_mutex_lock(&mutex);

  for (list<Hypothesis>::iterator it = hypotheses.begin(); it != hypotheses.end(); /*it++*/)
  {
    if (!is_p_value_significant(it->raw_log_p_value))
    {
      it->clustering_status = STATUS_REJECTED;
      it->M1_paired_matches.clear();
    }

    if (!is_p_value_significant(it->removal_raw_log_p_value))
    {
      hypotheses.erase(it++);
      hypotheses_size--;
    }
    else
      it++;
  }

  hypotheses_previous_size = hypotheses_size;
  pthread_mutex_unlock(&mutex);
  pthread_mutex_unlock(&remove_mutex);
}

long int HypothesesSet::getNumberOfAcceptedHypotheses() const
{
  long int count = 0;

  for (list<Hypothesis>::const_iterator it = hypotheses.begin(); it != hypotheses.end(); it++)
    if (it->clustering_status != STATUS_REJECTED)
      count++;

  return count;
}

void HypothesesSet::joinHypothesis(const Hypothesis *pred, Hypothesis *it, int offset, bool same_orientation, int clustering_status)
{
  it->cluster_id = pred->cluster_id;

  if (pred->cluster_same_orientation)
    it->cluster_offset = pred->cluster_offset + offset;
  else
    it->cluster_offset = pred->cluster_offset + pred->pair_end - pred->pair_start - it->pair_end + it->pair_start - offset;

  it->cluster_same_orientation = same_orientation ? pred->cluster_same_orientation : !pred->cluster_same_orientation;
  it->clustering_status = clustering_status;
}

void HypothesesSet::clusterHypotheses()
{
  for (list<Hypothesis>::iterator it = hypotheses.begin(); it != hypotheses.end(); it++)
    if (it->clustering_status != STATUS_REJECTED)
      clustered_hypotheses.push_back(&*it);

  sort(clustered_hypotheses.begin(), clustered_hypotheses.end(), Hypothesis::ptr_compare);

  for (vector<Hypothesis *>::iterator it_ptr = clustered_hypotheses.begin(); it_ptr != clustered_hypotheses.end(); it_ptr++)
  {
    Hypothesis *it = *it_ptr;

    if (spec.Options.ClusteringByIdentity || spec.Options.ClusteringDistanceConstant > 0.0 || spec.Options.ClusteringDistanceMultiplier > 0.0 || spec.Options.ClusteringOverlapThreshold < 1.0)
    {
      // try to join the overrepresented motif complex to any of the previous clusters
      for (vector<Hypothesis *>::iterator pred_ptr = clustered_hypotheses.begin(); pred_ptr != it_ptr; pred_ptr++)
      {
        Hypothesis *pred = *pred_ptr;

        if (!spec.Options.ClusteringAcrossDatasets && pred->dataset_id != it->dataset_id) continue;

        // first step of clustering: joining by motif complex identity
        // join only with the most significant motif complex in a previously established cluster
        if (spec.Options.ClusteringByIdentity)
        {
          if (pred->clustering_status == STATUS_CLUSTER_SEED)
          {
            if (pred->M1 == it->M1 && pred->M2 == it->M2 && pred->offset == it->offset && pred->same_orientation == it->same_orientation)
            {
              joinHypothesis(pred, it, 0, true, STATUS_JOINED_BY_IDENTITY);
              break;
            }
          }
        }

        // second step of clustering: joining by dimer motif similarity
        // join only with the signature motif complexes, i.e. identical to the most significant one in a previously established cluster
        if (spec.Options.ClusteringDistanceConstant > 0.0 || spec.Options.ClusteringDistanceMultiplier > 0.0)
        {
          if (pred->clustering_status == STATUS_CLUSTER_SEED || pred->clustering_status == STATUS_JOINED_BY_IDENTITY)
          {
            int offset;
            bool same_orientation;
            double dist = pred->dimer.getMinimumEuclideanDistanceSquared(it->dimer, &offset, &same_orientation);

            if (dist < spec.Options.ClusteringDistanceConstant + pred->dimer.inf_content * spec.Options.ClusteringDistanceMultiplier)
            {
              joinHypothesis(pred, it, offset, same_orientation, STATUS_JOINED_BY_SIMILARITY);
              break;
            }
          }
        }

        // third step of clustering: joining by overlap of genomic matches
        // (precisely, the overlap is calculated for the most frequent structure of the pair of motif dimers)
        // join only with the signature motif complexes and the complexes called similar to them at the second level
        if (spec.Options.ClusteringOverlapThreshold < 1.0)
        {
          if (pred->clustering_status == STATUS_CLUSTER_SEED || pred->clustering_status == STATUS_JOINED_BY_IDENTITY || pred->clustering_status == STATUS_JOINED_BY_SIMILARITY)
          {
            StatsSet ss(&*pred, &*it, spec.Options.MaxMotifSpacing);
            ss.calculateStats(pred->M1_paired_matches, it->M1_paired_matches);

            int offset;
            bool same_orientation;
            long int hits = ss.getMaximumStats(&offset, &same_orientation);

            if (hits >= spec.Options.ClusteringOverlapThreshold * max((double) it->target_instances - it->prob * it->target_N, 0.0))
            {
              joinHypothesis(pred, it, offset, same_orientation, STATUS_JOINED_BY_OVERLAP);
              break;
            }
          }
        }
      }
    }

    // if not joined, establish a new cluster
    if (it->clustering_status == STATUS_UNRESOLVED)
    {
      it->cluster_id = ++num_clusters;
      it->cluster_offset = 0;
      it->cluster_same_orientation = true;
      it->clustering_status = STATUS_CLUSTER_SEED;
    }
  }

  if (spec.Options.ClusteringByIdentity || spec.Options.ClusteringDistanceConstant > 0.0 || spec.Options.ClusteringDistanceMultiplier > 0.0 || spec.Options.ClusteringOverlapThreshold < 1.0)
    sort(clustered_hypotheses.begin(), clustered_hypotheses.end(), Hypothesis::ptr_compare);
}

void HypothesesSet::writeDetailedStatsFile(const char *fname) const
{
  if (spec.Options.OutputDetailedStats == OUTPUT_RANGE_NONE) return;
  cout << "Writing detailed stats file: " << fname << endl;
  FILE *fout;

  if ((fout = fopen(fname, "w")) == NULL)
  {
    cerr << "Could not open output file \"" << fname << "\"" << endl;
    exit(1);
  }

  fprintf(fout, "hypothesis_id\tM1_acc\tM2_acc\tdataset\toffset\torientation\ttarget_instances\ttarget_N\tcontrol_instances\tcontrol_N\tlog10_prob\tfold_change\tM1_inf_content\tM2_inf_content\toverlap_inf_content\tM1_inf_contribution\tM2_inf_contribution\tlog10_raw_p_value\tlog10_p_value\n");

  for (list<Hypothesis>::const_iterator it = hypotheses.begin(); it != hypotheses.end(); it++)
  {
    if (spec.Options.OutputDetailedStats == OUTPUT_RANGE_SIGNATURE && it->removal_hypothesis->clustering_status != STATUS_CLUSTER_SEED && it->removal_hypothesis->clustering_status != STATUS_JOINED_BY_IDENTITY) continue;

    fprintf(fout, "%ld\t%s\t%s\t%s\t%d\t%s\t%ld\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
      it->hypothesis_id, it->M1->accession.c_str(), it->M2->accession.c_str(),
      Genome::dataset_names[it->dataset_id].c_str(), it->offset, (it->same_orientation ? "same" : "opposite"),
      it->target_instances, it->target_N, it->control_instances, it->control_N, log(it->prob) / log(10), it->fold_change,
      it->M1->inf_content, it->M2->inf_content, it->overlap_inf_content, it->M1_inf_contribution, it->M2_inf_contribution,
      get_log10_raw_p_value(*it), get_log10_p_value(*it));
  }

  fclose(fout);
}

void HypothesesSet::writeGenomicLocationsFile(const char *fname) const
{
  if (spec.Options.OutputGenomicLocations == OUTPUT_RANGE_NONE) return;
  cout << "Writing genomic locations file: " << fname << endl;
  FILE *fout;

  if ((fout = fopen(fname, "w")) == NULL)
  {
    cerr << "Could not open output file \"" << fname << "\"" << endl;
    exit(1);
  }

  fprintf(fout, "hypothesis_id%s\tchrom\tstart\tend\tstrand\n", (spec.Options.GenomicLocationsMaxSpacingDeviation > 0) ? "\tspacing_deviation" : "");

  for (list<Hypothesis>::const_iterator it = hypotheses.begin(); it != hypotheses.end(); it++)
  {
    if (spec.Options.OutputGenomicLocations == OUTPUT_RANGE_SIGNATURE && it->clustering_status != STATUS_CLUSTER_SEED && it->clustering_status != STATUS_JOINED_BY_IDENTITY) continue;

    for (vector<PositionWeightMatrix::MotifMatch>::const_iterator jt = it->M1_paired_matches.begin(); jt != it->M1_paired_matches.end(); jt++)
      fprintf(fout, "%ld%s\t%s\t%d\t%d\t%c\n", it->hypothesis_id, (spec.Options.GenomicLocationsMaxSpacingDeviation > 0) ? "\t0" : "",
        Genome::chrom_names[jt->chrom_id].c_str(), it->get_start(*jt), it->get_end(*jt), jt->strand);

    int spacing_sign = it->M1_goes_first ? 1 : -1;
    for (int spacing_deviation = 1; spacing_deviation <= spec.Options.GenomicLocationsMaxSpacingDeviation; spacing_deviation++)
    {
      int spaced_pair_start = min(0, it->offset + spacing_sign * spacing_deviation);
      int spaced_pair_end = max(it->M1->length, it->M2->length + it->offset + spacing_sign * spacing_deviation);

      map<int, vector<PositionWeightMatrix::MotifMatch> >::const_iterator jt = it->M1_spaced_matches.find(spacing_deviation);
      if (jt != it->M1_spaced_matches.end())
        for (vector<PositionWeightMatrix::MotifMatch>::const_iterator kt = jt->second.begin(); kt != jt->second.end(); kt++)
          fprintf(fout, "%ld\t%d\t%s\t%d\t%d\t%c\n", it->hypothesis_id, spacing_deviation, Genome::chrom_names[kt->chrom_id].c_str(),
            (kt->strand == '+') ? kt->start + spaced_pair_start : kt->start + it->M1->length - spaced_pair_end,
            (kt->strand == '+') ? kt->start + spaced_pair_end : kt->start + it->M1->length - spaced_pair_start, kt->strand);
    }
  }

  fclose(fout);
}

void HypothesesSet::writePValueDistributionFile(const char *fname) const
{
  if (!spec.Options.OutputPValueDistribution) return;
  cout << "Writing p-value distribution file: " << fname << endl;
  FILE *fout;

  if ((fout = fopen(fname, "w")) == NULL)
  {
    cerr << "Could not open output file \"" << fname << "\"" << endl;
    exit(1);
  }

  fprintf(fout, "log10_raw_p_value\tcount\n");

  for (int i = 0; i < (int) log10_raw_p_value_dist.size(); i++)
    fprintf(fout, "%f\t%ld\n", -(i + 1) * log10_raw_p_value_sampling, log10_raw_p_value_dist[i]);

  fclose(fout);
}

void HypothesesSet::writeClusteringResultsFile(const char *fname) const
{
  FILE *fout;
  cout << "Writing clustering results file: " << fname << endl;

  if ((fout = fopen(fname, "w")) == NULL)
  {
    cerr << "Could not open output file \"" << fname << "\"" << endl;
    exit(1);
  }

  fprintf(fout, "cluster_id\thypothesis_id\tM1_acc\tM1_name\tM1_orientation\toffset\tM2_acc\tM2_name\tM2_orientation\tdataset\ttarget_instances\ttarget_N\tcontrol_instances\tcontrol_N\tlog10_prob\tfold_change\tM1_inf_content\tM2_inf_content\toverlap_inf_content\tM1_inf_contribution\tM2_inf_contribution\tlog10_raw_p_value\tlog10_p_value\tcluster_offset\tsimilarity_annotation\n");

  for (vector<Hypothesis *>::const_iterator it_ptr = clustered_hypotheses.begin(); it_ptr != clustered_hypotheses.end(); it_ptr++)
  {
    Hypothesis *it = *it_ptr;

    if (it->clustering_status == STATUS_REJECTED) continue;

    bool M1_goes_first, M1_same_orientation, M2_same_orientation, similarity_same_orientation;
    int offset, similarity_offset;

    if (it->cluster_same_orientation)
    {
      M1_goes_first = it->M1_goes_first;
      M1_same_orientation = true;
      M2_same_orientation = it->same_orientation;
      offset = it->offset;
      similarity_same_orientation = it->similarity.same_orientation;
      similarity_offset = it->similarity.offset;
    }
    else
    {
      M1_goes_first = !it->M1_goes_first;
      M1_same_orientation = false;
      M2_same_orientation = !it->same_orientation;
      offset = it->M1->length - it->offset - it->M2->length;
      similarity_same_orientation = !it->similarity.same_orientation;
      similarity_offset = it->pair_end - it->pair_start - it->similarity.offset - it->similarity.motif->length;
    }

    fprintf(fout, "%d\t%ld\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%ld\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t",
      it->cluster_id, it->hypothesis_id,
      M1_goes_first ? it->M1->accession.c_str() : it->M2->accession.c_str(),
      M1_goes_first ? it->M1->name.c_str() : it->M2->name.c_str(),
      M1_goes_first ? ots(M1_same_orientation) : ots(M2_same_orientation),
      M1_goes_first ? offset : -offset,
      M1_goes_first ? it->M2->accession.c_str() : it->M1->accession.c_str(),
      M1_goes_first ? it->M2->name.c_str() : it->M1->name.c_str(),
      M1_goes_first ? ots(M2_same_orientation) : ots(M1_same_orientation),
      Genome::dataset_names[it->dataset_id].c_str(),
      it->target_instances, it->target_N, it->control_instances, it->control_N, log(it->prob) / log(10), it->fold_change,
      M1_goes_first ? it->M1->inf_content : it->M2->inf_content,
      M1_goes_first ? it->M2->inf_content : it->M1->inf_content,
      it->overlap_inf_content,
      M1_goes_first ? it->M1_inf_contribution : it->M2_inf_contribution,
      M1_goes_first ? it->M2_inf_contribution : it->M1_inf_contribution,
      get_log10_raw_p_value(*it), get_log10_p_value(*it), it->cluster_offset);

    if (it->similarity.distance < spec.Options.ClusteringDistanceConstant + it->dimer.inf_content * spec.Options.ClusteringDistanceMultiplier)
      fprintf(fout, "<- similar to %s (%s) [%d, %s, ED2 %g]", it->similarity.motif->accession.c_str(), it->similarity.motif->name.c_str(), similarity_offset, (similarity_same_orientation ? "same" : "opposite"), it->similarity.distance);

    fprintf(fout, "\n");
  }

  fclose(fout);
}

void HypothesesSet::writeDimerMotifsFile(const char *fname) const
{
  if (spec.Options.OutputDimerMotifs == OUTPUT_RANGE_NONE) return;
  cout << "Writing dimer motifs file: " << fname << endl;
  FILE *fout;

  if ((fout = fopen(fname, "w")) == NULL)
  {
    cerr << "Could not open output file \"" << fname << "\"" << endl;
    exit(1);
  }

  int cluster_seq = 0;
  for (vector<Hypothesis *>::const_iterator it_ptr = clustered_hypotheses.begin(); it_ptr != clustered_hypotheses.end(); it_ptr++)
  {
    Hypothesis *it = *it_ptr;

    if (it->clustering_status == STATUS_REJECTED) continue;
    if (it->clustering_status == STATUS_CLUSTER_SEED) cluster_seq = 0;
    cluster_seq++;

    if (spec.Options.OutputDimerMotifs == OUTPUT_RANGE_SIGNATURE && it->clustering_status != STATUS_CLUSTER_SEED && it->clustering_status != STATUS_JOINED_BY_IDENTITY) continue;

    bool M1_goes_first, M1_same_orientation, M2_same_orientation;
    int offset;

    if (it->cluster_same_orientation)
    {
      M1_goes_first = it->M1_goes_first;
      M1_same_orientation = true;
      M2_same_orientation = it->same_orientation;
      offset = it->offset;
    }
    else
    {
      M1_goes_first = !it->M1_goes_first;
      M1_same_orientation = false;
      M2_same_orientation = !it->same_orientation;
      offset = it->M1->length - it->offset - it->M2->length;
    }

    fprintf(fout, "AC  cluster_%d_motif_complex_%d\nXX\nID  %s_%s_%d_%s_%s_%s\nXX\nNA  hypothesis_%ld\nXX\nP0\tA\tC\tG\tT\n",
      it->cluster_id, cluster_seq,
      M1_goes_first ? it->M1->accession.c_str() : it->M2->accession.c_str(),
      M1_goes_first ? ots(M1_same_orientation) : ots(M2_same_orientation),
      M1_goes_first ? offset : -offset,
      M1_goes_first ? it->M2->accession.c_str() : it->M1->accession.c_str(),
      M1_goes_first ? ots(M2_same_orientation) : ots(M1_same_orientation),
      Genome::dataset_names[it->dataset_id].c_str(),
      it->hypothesis_id);

    if (it->cluster_same_orientation)
      for (int i = 0; i < it->dimer.length; i++)
        fprintf(fout, "%02d\t%f\t%f\t%f\t%f\n", i,
          it->dimer.probabilities[i].values[nti('A')], it->dimer.probabilities[i].values[nti('C')],
          it->dimer.probabilities[i].values[nti('G')], it->dimer.probabilities[i].values[nti('T')]);
    else
      for (int i = it->dimer.length - 1; i >= 0; i--)
        fprintf(fout, "%02d\t%f\t%f\t%f\t%f\n", it->dimer.length - i - 1,
          it->dimer.probabilities[i].values[nti('T')], it->dimer.probabilities[i].values[nti('G')],
          it->dimer.probabilities[i].values[nti('C')], it->dimer.probabilities[i].values[nti('A')]);

    fprintf(fout, "XX\n//\n");
  }

  fclose(fout);
}

HypothesesSet::~HypothesesSet()
{
  pthread_mutex_destroy(&mutex);
  pthread_mutex_destroy(&remove_mutex);
}
