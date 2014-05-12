/*
    PositionWeightMatrix.cpp

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
#include <math.h>
#include <sstream>
#include <string.h>

#include "PositionWeightMatrix.h"

#define FORMAT_UNKNOWN 0
#define FORMAT_TRANSFAC 1
#define FORMAT_JASPAR 2
#define FORMAT_MEME 3

const int score_samples = 1000;
const double pseudocounts = 0.05;

#define nti(nucleotide) Genome::nucleotideToIndex(nucleotide)

static inline double square(double x)
{
  return x * x;
}

inline void PositionWeightMatrix::initialize(double GC_content)
{
  length = 0;
  inf_content = 0.0;
  bgprob.values[0] = 0.0;
  bgprob.values[nti('A')] = (1 - GC_content) / 2;
  bgprob.values[nti('C')] = GC_content / 2;
  bgprob.values[nti('G')] = GC_content / 2;
  bgprob.values[nti('T')] = (1 - GC_content) / 2;
  threshold = -1;
}

PositionWeightMatrix::PositionWeightMatrix(double GC_content)
{
  initialize(GC_content);
}

PositionWeightMatrix::PositionWeightMatrix(double GC_content, const char *fname, const string &accession)
{
  initialize(GC_content);
  if (accession != "") name = accession;

  FILE *f = fopen(fname, "r");
  if (f == NULL)
  {
    cerr << fname << ": " << strerror(errno) << endl;
    exit(1);
  }
  int lineNum = 0;

  try
  {
    readPWM(f, fname, &lineNum);
  }
  catch (PositionWeightMatrix::EndOfFile &e)
  {
    cerr << fname << ": premature end of file" << endl;
    exit(1);
  }

  fclose(f);

  if (accession != "") this->accession = accession;
}

PositionWeightMatrix::PositionWeightMatrix(double GC_content, FILE *f, const char *fname, int *lineNum)
{
  initialize(GC_content);
  readPWM(f, fname, lineNum);
}

void PositionWeightMatrix::readPWM(FILE *f, const char *fname, int *lineNum)
{
  PWM.clear();
  char line[1024];
  int format = FORMAT_UNKNOWN;
  bool first_line = false, readingPWM = false;

  // guess the PWM format: TRANSFAC or JASPAR?
  while (fgets(line, 1024, f))
  {
    (*lineNum)++;
    istringstream is(line);
    string label;
    if (is >> label)
    {
      if (label[0] == '>')
        format = FORMAT_JASPAR;
      else if (label == "letter-probability")
        format = FORMAT_MEME;
      else
        format = FORMAT_TRANSFAC;

      first_line = true;
      break;
    }
  }

  if (format == FORMAT_TRANSFAC) // also SwissRegulon
  {
    int index[4] = {0, 0, 0, 0};

    while (first_line || fgets(line, 1024, f))
    {
      if (first_line) first_line = false; else (*lineNum)++;
      istringstream is(line);
      string label;
      is >> label >> ws;

      if (label == "AC") getline(is, accession);
      if (label == "ID") getline(is, id);
      if (label == "NA") getline(is, name);

      if (label == "P0" || label == "PO" || label == "POS")
      {
        readingPWM = true;
        for (int i = 0; i < 4; i++)
        {
          char nucleotide;
          if (!(is >> nucleotide))
          {
            cerr << fname << ": not enough nucleotide labels in line " << *lineNum << endl;
            exit(1);
          }
          index[i] = nti(nucleotide);
        }
        continue;
      }

      if (readingPWM)
      {
        if (label == "XX" || label == "//") break;

        NucleotideValues fd;
        for (int i = 0; i < 5; i++)
          fd.values[i] = 0.0;

        for (int i = 0; i < 4; i++)
        {
          if (!(is >> fd.values[index[i]]))
          {
            cerr << fname << ": not enough values in line " << *lineNum << endl;
            exit(1);
          }
        }

        if (addPosition(fd, pseudocounts) != 0)
        {
          cerr << fname << ": malformed nucleotide weights in line " << *lineNum << endl;
          exit(1);
        }
      }
    }
  }

  if (format == FORMAT_JASPAR)
  {
    istringstream is(line);
    char c;
    is >> c >> accession >> name;

    first_line = false;
    readingPWM = true;
    map<int, vector<double> > values;

    int gc;
    while ((gc = getc(f)) != EOF)
    {
      ungetc(gc, f);
      if (gc == '>') break;
      if (!fgets(line, 1024, f)) break;
      (*lineNum)++;

      istringstream is(line);
      char nucleotide;
      if (!(is >> nucleotide)) continue;
      int index = nti(nucleotide);

      if (!values[index].empty())
      {
        cerr << fname << ": repeat of nucleotide " << nucleotide << " in line " << *lineNum << endl;
        exit(1);
      }

      is >> ws;
      if (is.peek() == '[') is.get();
      double value;
      while (is >> value)
        values[index].push_back(value);
    }

    if (values.size() == 4)
    {
      int my_length = values.begin()->second.size();
      for (map<int, vector<double> >::iterator it = values.begin(); it != values.end(); it++)
        if ((int) it->second.size() != my_length)
        {
          cerr << fname << ": rows of different length in line " << *lineNum << endl;
          exit(1);
        }

      for (int i = 0; i < my_length; i++)
      {
        NucleotideValues fd;
        for (int j = 0; j < 5; j++)
          fd.values[j] = 0.0;

        for (map<int, vector<double> >::iterator it = values.begin(); it != values.end(); it++)
          fd.values[it->first] = it->second[i];

        if (addPosition(fd, pseudocounts) != 0)
        {
          cerr << fname << ": malformed nucleotide weights in line " << *lineNum << endl;
          exit(1);
        }
      }
    }
  }

  if (format == FORMAT_MEME)
  {
    while (fgets(line, 1024, f))
    {
      first_line = false;
      readingPWM = true;
      (*lineNum)++;

      istringstream is(line);
      NucleotideValues fd;
      for (int i = 0; i < 5; i++)
        fd.values[i] = 0.0;

      for (int i = 1; i < 5; i++)
      {
        if (!(is >> fd.values[i]))
        {
          if (i > 0)
          {
            cerr << fname << ": not enough values in line " << *lineNum << endl;
            exit(1);
          }
          readingPWM = false;
          break;
        }
      }

      if (readingPWM)
      {
        if (addPosition(fd, pseudocounts) != 0)
        {
          cerr << fname << ": malformed nucleotide weights in line " << *lineNum << endl;
          exit(1);
        }
      }
      else break;
    }
  }

  if (format == FORMAT_UNKNOWN && !feof(f))
  {
    cerr << fname << ": unknown PWM format in line " << *lineNum << endl;
    exit(1);
  }

  if (length == 0 && feof(f)) throw EndOfFile();

  if (accession == "")
  {
    if (id != "") accession = id;
    else if (name != "") accession = name;
    else accession = fname;
  }
}

int PositionWeightMatrix::addPosition(NucleotideValues fd, double pseudocounts)
{
  length++;

  double sum = 0.0;
  for (int i = 1; i < 5; i++)
  {
    if (fd.values[i] < 0) return E_NEGATIVE_WEIGHT;
    sum += fd.values[i];
  }
  if (sum == 0) return E_ALL_WEIGHTS_ZERO;

  // normalize the frequency distribution, add pseudo-counts
  fd.values[nti('A')] = (1 - pseudocounts) * fd.values[nti('A')] / sum + pseudocounts * bgprob.values[nti('A')];
  fd.values[nti('C')] = (1 - pseudocounts) * fd.values[nti('C')] / sum + pseudocounts * bgprob.values[nti('C')];
  fd.values[nti('G')] = (1 - pseudocounts) * fd.values[nti('G')] / sum + pseudocounts * bgprob.values[nti('G')];
  fd.values[nti('T')] = (1 - pseudocounts) * fd.values[nti('T')] / sum + pseudocounts * bgprob.values[nti('T')];
  probabilities.push_back(fd);

  // update the information content
  double fd_ic = 0.0;
  if (fd.values[nti('A')] > 0.0) fd_ic += fd.values[nti('A')] * (log(fd.values[nti('A')]) - log(bgprob.values[nti('A')])) / log(2);
  if (fd.values[nti('C')] > 0.0) fd_ic += fd.values[nti('C')] * (log(fd.values[nti('C')]) - log(bgprob.values[nti('C')])) / log(2);
  if (fd.values[nti('G')] > 0.0) fd_ic += fd.values[nti('G')] * (log(fd.values[nti('G')]) - log(bgprob.values[nti('G')])) / log(2);
  if (fd.values[nti('T')] > 0.0) fd_ic += fd.values[nti('T')] * (log(fd.values[nti('T')]) - log(bgprob.values[nti('T')])) / log(2);
  position_inf_content.push_back(fd_ic);
  inf_content += fd_ic;

  // calculate log-odds wrt to the background distribution
  fd.values[nti('A')] = (log(fd.values[nti('A')]) - log(bgprob.values[nti('A')])) / log(2);
  fd.values[nti('C')] = (log(fd.values[nti('C')]) - log(bgprob.values[nti('C')])) / log(2);
  fd.values[nti('G')] = (log(fd.values[nti('G')]) - log(bgprob.values[nti('G')])) / log(2);
  fd.values[nti('T')] = (log(fd.values[nti('T')]) - log(bgprob.values[nti('T')])) / log(2);
  fd.values[0] = fd.values[nti('A')] * bgprob.values[nti('A')] + fd.values[nti('C')] * bgprob.values[nti('C')]
    + fd.values[nti('G')] * bgprob.values[nti('G')] + fd.values[nti('T')] * bgprob.values[nti('T')];
  PWM.push_back(fd);

  return 0;
}

void PositionWeightMatrix::estimateThreshold(double sensitivity)
{
  min_score = max_score = 0.0;

  for (int i = 0; i < length; i++)
  {
    min_score += min(min(PWM[i].values[nti('A')], PWM[i].values[nti('C')]), min(PWM[i].values[nti('G')], PWM[i].values[nti('T')]));
    max_score += max(max(PWM[i].values[nti('A')], PWM[i].values[nti('C')]), max(PWM[i].values[nti('G')], PWM[i].values[nti('T')]));
  }

  double score_probs_carrier[score_samples + 1];
  double *score_probs = score_probs_carrier;
  double score_probs_new_carrier[score_samples + 1];
  double *score_probs_new = score_probs_new_carrier;

  for (int i = 0; i < score_samples + 1; i++)
    score_probs[i] = 0.0;
  score_probs[(int) round((0 - min_score) / (max_score - min_score) * score_samples)] = 1.0;

  for (int i = 0; i < length; i++)
  {
    for (int j = 0; j < score_samples + 1; j++)
      score_probs_new[j] = 0.0;
    for (int k = 1; k < 5; k++)
    {
      int score_incr = (int) round(PWM[i].values[k] / (max_score - min_score) * score_samples);
      for (int j = 0; j < score_samples + 1; j++)
        score_probs_new[min(max(j + score_incr, 0), score_samples)] += score_probs[j] * probabilities[i].values[k];
    }
    swap(score_probs, score_probs_new);
  }

  double sum = 0;
  int index = score_samples;
  while (sum < sensitivity && index > 0)
    sum += score_probs[index--];
  threshold = min_score + index * (max_score - min_score) / score_samples;

  // now assuming background distribution
  for (int i = 0; i < score_samples + 1; i++)
    score_probs[i] = 0.0;
  score_probs[(int) round((0 - min_score) / (max_score - min_score) * score_samples)] = 1.0;

  for (int i = 0; i < length; i++)
  {
    for (int j = 0; j < score_samples + 1; j++)
      score_probs_new[j] = 0.0;
    for (int k = 1; k < 5; k++)
    {
      int score_incr = (int) round(PWM[i].values[k] / (max_score - min_score) * score_samples);
      for (int j = 0; j < score_samples + 1; j++)
        score_probs_new[min(max(j + score_incr, 0), score_samples)] += score_probs[j] * bgprob.values[k];
    }
    swap(score_probs, score_probs_new);
  }

  sum = 0;
  for (int i = 0; i < index + 1; i++)
    sum += score_probs[i];
  cout << "Given sensitivity " << sensitivity << ", calculated threshold " << threshold << " and specificity " << sum << endl;
}

double PositionWeightMatrix::forwardScore(const char *seq) const
{
  double score = 0.0;

  for (int i = 0; i < length; i++)
    score += PWM[i].values[nti(seq[i])];

  return score;
}

double PositionWeightMatrix::reverseScore(const char *seq) const
{
  double score = 0.0;

  for (int i = 0; i < length; i++)
    score += PWM[i].values[Genome::complementIndex(nti(seq[length - 1 - i]))];

  return score;
}

void PositionWeightMatrix::scan(const Genome &fa, NarrowPeak *regions)
{
  matches.clear();

  if (regions)
    for (list<NarrowPeak::Region>::iterator jt = regions->regions.begin(); jt != regions->regions.end(); jt++)
    {
      map<int, vector<char> >::const_iterator it = fa.chroms.find(jt->chrom_id);
      if (it == fa.chroms.end()) continue;
      // pseudo_start and pseudo_end set the scope for leftmost PWM postions, cf. StatsSet::addRegions
      int jt_pseudo_start = max(jt->start - length / 2, 0);
      int jt_pseudo_end = min(jt->end - length / 2, (int) it->second.size() - length + 1);

      for (int j = jt_pseudo_start; j < jt_pseudo_end; j++)
      {
        if (forwardScore(&it->second[j]) > threshold)
        {
          MotifMatch mm;
          mm.chrom_id = it->first;
          mm.start = j;
          mm.strand = '+';
          matches.push_back(mm);
        }

        if (reverseScore(&it->second[j]) > threshold)
        {
          MotifMatch mm;
          mm.chrom_id = it->first;
          mm.start = j;
          mm.strand = '-';
          matches.push_back(mm);
        }
      }
    }
  else
    for (map<int, vector<char> >::const_iterator it = fa.chroms.begin(); it != fa.chroms.end(); it++)
    {
      for (unsigned int j = 0; j < it->second.size() - length + 1; j++)
      {
        if (forwardScore(&it->second[j]) > threshold)
        {
          MotifMatch mm;
          mm.chrom_id = it->first;
          mm.start = j;
          mm.strand = '+';
          matches.push_back(mm);
        }

        if (reverseScore(&it->second[j]) > threshold)
        {
          MotifMatch mm;
          mm.chrom_id = it->first;
          mm.start = j;
          mm.strand = '-';
          matches.push_back(mm);
        }
      }
    }

  // shrink the matches vector, since it will not be expanded any more
  vector<MotifMatch>(matches).swap(matches);
  // and sort it (by the default method, by chrom_id and start)
  sort(matches.begin(), matches.end());
}

bool PositionWeightMatrix::MotifMatch::operator<(const MotifMatch &other) const
{
  if (this->chrom_id < other.chrom_id) return true;
  else if (this->chrom_id > other.chrom_id) return false;
  else return this->start < other.start;
}

double PositionWeightMatrix::EuclideanDistanceSquared(const PositionWeightMatrix &other, int offset, bool same_orientation) const
{
  int pair_start = min(0, offset);
  int pair_end = max(this->length, other.length + offset);
  double result = 0.0;

  for (int i = pair_start; i < pair_end; i++)
  {
    const NucleotideValues *p, *q;

    if (i >= 0 && i < this->length)
      p = &this->probabilities[i];
    else
      p = &this->bgprob;

    if (i >= offset && i < other.length + offset)
      if (same_orientation)
        q = &other.probabilities[i - offset];
      else
        q = &other.probabilities[other.length - 1 - i + offset];
    else
      q = &other.bgprob;

    if (same_orientation)
      result += square(p->values[nti('A')] - q->values[nti('A')]) + square(p->values[nti('C')] - q->values[nti('C')])
        + square(p->values[nti('G')] - q->values[nti('G')]) + square(p->values[nti('T')] - q->values[nti('T')]);
    else
      result += square(p->values[nti('A')] - q->values[nti('T')]) + square(p->values[nti('C')] - q->values[nti('G')])
        + square(p->values[nti('G')] - q->values[nti('C')]) + square(p->values[nti('T')] - q->values[nti('A')]);
  }

  return result;
}

double PositionWeightMatrix::getMinimumEuclideanDistanceSquared(const PositionWeightMatrix &other, int *returned_offset, bool *returned_same_orientation) const
{
  double min_dist = INFINITY;

  int offset;
  bool same_orientation = true;
  while (true) // iterating same_orientation over {true, false}
  {
    for (offset = -other.length + 1; offset < this->length; offset++)
    {
      double dist = this->EuclideanDistanceSquared(other, offset, same_orientation);

      if (dist < min_dist)
      {
        min_dist = dist;
        *returned_offset = offset;
        *returned_same_orientation = same_orientation;
      }
    }

    if (same_orientation) same_orientation = false; else break;
  }

  return min_dist;
}
