/*
    HypothesesSet.h

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

#ifndef _HYPOTHESES_SET_H
#define _HYPOTHESES_SET_H

#include <list>
#include <math.h>

#include "PositionWeightMatrix.h"
#include "Specification.h"

using namespace std;

// initial values for hypothesis_id
#define HYPOTHESIS_CONSIDERED -2
#define HYPOTHESIS_NOT_CONSIDERED -1

// clustering status
#define STATUS_UNRESOLVED 0
#define STATUS_REJECTED -1
#define STATUS_JOINED_BY_IDENTITY 1
#define STATUS_JOINED_BY_SIMILARITY 2
#define STATUS_JOINED_BY_OVERLAP 3
#define STATUS_CLUSTER_SEED 4

class HypothesesSet
{
  public:
    struct Hypothesis
    {
      long int hypothesis_id, target_instances, target_N, control_instances, control_N;
      double freq_ratio, prob, fold_change, overlap_inf_content, M1_inf_contribution, M2_inf_contribution;
      double raw_log_p_value, removal_raw_log_p_value;

      PositionWeightMatrix *M1, *M2;
      PositionWeightMatrix dimer;
      Hypothesis *removal_hypothesis;

      vector<PositionWeightMatrix::MotifMatch> M1_paired_matches;
      map<int, vector<PositionWeightMatrix::MotifMatch> > M1_spaced_matches;

      struct
      {
        double distance;
        const PositionWeightMatrix *motif;
        int offset;
        bool same_orientation;
      } similarity;

      int dataset_id, offset, pair_start, pair_end;
      int clustering_status, cluster_id, cluster_offset;
      bool same_orientation, M1_goes_first, cluster_same_orientation;

      bool operator<(const Hypothesis &other) const;
      static bool ptr_compare(const Hypothesis *lhs, const Hypothesis *rhs);

      void calculateStructure();
      void calculateDimerMotif(const Genome &fa, double GC_content, int dimer_motif_margin);
      void calculateDimerMotifSimilarity(const Genome &fa, const vector<PositionWeightMatrix> &PWMs, double GC_content, int dimer_motif_margin = 0);
      int get_start(const PositionWeightMatrix::MotifMatch &mm) const; // start and end of a particular motif pair instance
      int get_end(const PositionWeightMatrix::MotifMatch &mm) const;
      int get_length() const; // length in bp
    };

    const Specification &spec;
    HypothesesSet(const Specification &spec);
    void processPair(PositionWeightMatrix *M1, PositionWeightMatrix *M2, const map<int, const NarrowPeak *> &target_datasets, const NarrowPeak *control_dataset);
    void removeInsignificantHypotheses();
    long int getNumberOfHypotheses() const;
    long int getNumberOfAcceptedHypotheses() const;
    void clusterHypotheses();
    int getNumberOfClusters() const;

    void writeDetailedStatsFile(const char *fname) const;
    void writeGenomicLocationsFile(const char *fname) const;
    void writePValueDistributionFile(const char *fname) const;
    void writeClusteringResultsFile(const char *fname) const;
    void writeDimerMotifsFile(const char *fname) const;

    double get_log10_raw_p_value(const Hypothesis &hyp) const;
    double get_log10_p_value(const Hypothesis &hyp) const;
    bool is_p_value_significant(double p_value) const;
    list<Hypothesis>::iterator hypotheses_begin();
    list<Hypothesis>::iterator hypotheses_end();
    ~HypothesesSet();

  private:
    double max_log_p_value;
    vector<long int> log10_raw_p_value_dist;

    long int num_hypotheses, hypotheses_size, hypotheses_previous_size;
    pthread_mutex_t data_mutex, remove_mutex;
    list<Hypothesis> hypotheses;

    int num_clusters;
    vector<Hypothesis *> clustered_hypotheses;
    void joinHypothesis(const Hypothesis *pred, Hypothesis *it, int offset, bool same_orientation, int clustering_status);
};

inline int HypothesesSet::Hypothesis::get_start(const PositionWeightMatrix::MotifMatch &mm) const
{
  return (mm.strand == '+') ? mm.start + pair_start : mm.start + M1->length - pair_end;
}

inline int HypothesesSet::Hypothesis::get_end(const PositionWeightMatrix::MotifMatch &mm) const
{
  return (mm.strand == '+') ? mm.start + pair_end : mm.start + M1->length - pair_start;
}

inline int HypothesesSet::Hypothesis::get_length() const
{
  return pair_end - pair_start;
}

inline long int HypothesesSet::getNumberOfHypotheses() const
{
  return num_hypotheses;
}

inline int HypothesesSet::getNumberOfClusters() const
{
  return num_clusters;
}

inline double HypothesesSet::get_log10_raw_p_value(const Hypothesis &h) const
{
  return h.raw_log_p_value / log(10);
}

inline double HypothesesSet::get_log10_p_value(const Hypothesis &h) const
{
  return (h.raw_log_p_value < INFINITY) ? min(h.raw_log_p_value + log(num_hypotheses), 0.0) / log(10) : INFINITY;
}

inline bool HypothesesSet::is_p_value_significant(double p_value) const
{
  return p_value + log(num_hypotheses) < max_log_p_value || max_log_p_value == INFINITY;
}

inline list<HypothesesSet::Hypothesis>::iterator HypothesesSet::hypotheses_begin()
{
  return hypotheses.begin();
}

inline list<HypothesesSet::Hypothesis>::iterator HypothesesSet::hypotheses_end()
{
  return hypotheses.end();
}

#endif
