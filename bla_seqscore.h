/*
  Copyright (c) 2012-2013 Meike Bruns, Florian Markowsky, Michael Spohn
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef BLA_SEQSCORE_H
#define BLA_SEQSCORE_H
#include "core/encseq.h"
#include "core/error_api.h"
#include "match/sarr-def.h"
#include "extended/ranked_list.h"

/* Struct that stores sequence numbers <seq_a> of the query and <seq_b> of the
   database sequence along with their <score>. GtRankedLists used below contain
   GtSeqPairs.*/
typedef struct
{
  unsigned long seq_a;
  unsigned long seq_b;
  double score;
} GtSeqPair;

/* Type definition of callback-function for KSEARCH results. The callback
   function gets sequence numbers <seqno1> and <seqno2> of query- and database-
   sequence along with their <score>, a void* <data> that may contain a result
   structure (if one is needed, has to be castet to its real data type) and an
   error object <err>.*/
typedef int (*ResultProcessor)(unsigned long seqno1, unsigned long seqno2,
                               unsigned long score, void *data, GtError *err);

/* Enumerates over KSEARCH scores for query <encseq_q> and database <encseq_s>.
   The wordsize for the k-mer profiles is <k>. For every sequence pair whose
   score is calculated, the callback function defined in <ResultProcessor> is
   called. Returns -1 on error, <err> is set accordingly.*/
int gt_seqscore_enumerate_results_ksearch(ResultProcessor, void *data,
                                          const GtEncseq* encseq_q,
                                          const GtEncseq* encseq_s,
                                          unsigned int k, GtError *err);
/* Performs SANS scoring for query <suf_q> and database <suf_s> with 
   windowsize <h>. Sequence pairs with a score higher than <threshold> are
   written to <filename> with their score.*/
void gt_seqscore_threshold_new_sans(Suffixarray suf_q, Suffixarray suf_s,
                                    unsigned int h, GtFile *filename,
                                    unsigned int threshold, GtError* err);

/* Returns a new <GtRankedList> containing those <nbest> pairs of sequence
   numbers from query <encseq_q> and database <encseq_s> with the highest
   KSEARCH scores. The wordsize for the k-mer profiles is <k>. Returns NULL on
   error, <err> is set accordingly. */
GtRankedList* gt_seqscore_ranked_list_new_ksearch(const GtEncseq *encseq_q,
                                                  const GtEncseq *encseq_s,
                                                  unsigned int k,
                                                  unsigned int nbest,
                                                  GtError *err);

/* Returns a new <GtRankedList> containing those <nbest> pairs of sequence
   numbers from query <suf_q> and database <suf_s> with the highest
   scores. The windowsize for the SANS algorithm is <h>. Returns NULL on
   error, <err> is set accordingly*/
GtRankedList* gt_seqscore_ranked_list_new_sans(Suffixarray suf_q,
                                               Suffixarray suf_s,
                                               unsigned int h,
                                               unsigned int nbest,
                                               GtError* err);

/* Definition of result matrix to store scores of every query/database sequence
   pair.*/
typedef double** GtSeqScoreResultMatrix;
/* Returns a <GtSeqScoreResultMatrix> containing all SANS scores of sequences
   of query <suf_q> and database <suf_s>. The windowsize for the SANS algorithm
   is <h>. Returns NULL on error, <err> is set accordingly.*/
GtSeqScoreResultMatrix gt_seqscore_result_matrix_new_sans(Suffixarray suf_q,
                                                          Suffixarray suf_s,
                                                          unsigned int h,
                                                          GtError* err);
/* Returns a <GtSeqScoreResultMatrix> containing all iKSEARCH scores of
   sequences of query <suf_q> and database <suf_s>. The wordsize for the k-mer
   profiles is <k>. Returns NULL on error, <err> is set accordingly.*/
GtSeqScoreResultMatrix gt_seqscore_result_matrix_new_ksearch(
                                                       const GtEncseq* encseq_q,
                                                       const GtEncseq* encseq_s,
                                                       unsigned int ksearch_k,
                                                       GtError* err);
/* Returns the score of query sequence number <seq_a> and database sequence
   number <seq_b> of <GtSeqScoreResultMatrix> <m>.*/
double                 gt_seqscore_result_matrix_get_score(
                                                       GtSeqScoreResultMatrix m,
                                                       unsigned long seq_a,
                                                       unsigned long seq_b);
/* Frees the allocated space of <GtSeqScoreResultMatrix> <m>.*/
void                   gt_seqscore_result_matrix_delete(GtSeqScoreResultMatrix);
#endif
