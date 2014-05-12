/*
    PositionWeightMatrix.h

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

#ifndef _POSITION_WEIGHT_MATRIX_H
#define _POSITION_WEIGHT_MATRIX_H

#include <string>
#include <vector>

#include "Genome.h"
#include "NarrowPeak.h"

using namespace std;

#define E_NEGATIVE_WEIGHT 1
#define E_ALL_WEIGHTS_ZERO 1

class PositionWeightMatrix
{
  public:
    string accession, id, name;
    int length;

    struct NucleotideValues
    {
      double values[5];
    };
    vector<NucleotideValues> probabilities, PWM;
    NucleotideValues bgprob;

    double min_score, max_score, threshold;
    double inf_content;
    vector<double> position_inf_content;

    struct MotifMatch
    {
      int start;
      short int chrom_id;
      char strand;
      bool operator<(const MotifMatch &other) const;
    };
    vector<MotifMatch> matches;

    PositionWeightMatrix(double GC_content = 0.25);
    PositionWeightMatrix(double GC_content, const char *fname, const string &accession = "");
    PositionWeightMatrix(double GC_content, FILE *f, const char *fname, int *lineNum);
    class EndOfFile : exception {};

    int addPosition(NucleotideValues fd, double pseudocounts = 0.0);
    void estimateThreshold(double sensitivity);
    double forwardScore(const char *seq) const;
    double reverseScore(const char *seq) const;
    void scan(const Genome &fa, NarrowPeak *regions = NULL);
    double EuclideanDistanceSquared(const PositionWeightMatrix &other, int offset, bool same_orientation) const;
    double getMinimumEuclideanDistanceSquared(const PositionWeightMatrix &other, int *returned_offset, bool *returned_same_orientation) const;

  protected:
    void initialize(double GC_content);
    void readPWM(FILE *f, const char *fname, int *lineNum);
};

#endif
