/*
 * Copyright [1999-2018] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
  Bamcount

  A program for calculating RPKM type values for features in a BAM format
  against a set of annotation from an Ensembl database

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
*/

#include <stdio.h>

#include "ensc/EnsC.h"

#include "ensc/DBAdaptor.h"
#include "ensc/BaseAdaptor.h"
#include "ensc/Vector.h"
#include "ensc/SliceAdaptor.h"
#include "ensc/Slice.h"
#include "ensc/DNAAlignFeature.h"
#include "ensc/StrUtil.h"
#include "ensc/IDHash.h"
#include "ensc/Transcript.h"

#include "bamhelper.h"
#include "htslib/sam.h"
#include "htslib/hts.h"

void       Bamcount_usage();
int        countReads(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags, Vector *genes, IDHash *geneResultsHash, long long countUsableReads, FILE *bedFp);

Vector *Bam_cigarToUngapped(bam1_t *b);

// Flag values
#define M_UCSC_NAMING 1
#define M_VERIFY 4

// My custom bam core struct flag value
#define MY_FUSEDFLAG 32768

int verbosity = 1;

int main(int argc, char *argv[]) {
  DBAdaptor *      dba;
  StatementHandle *sth;
  ResultRow *      row;
  Vector *         slices;
  int              nSlices;

  int   argNum = 1;

  char *inFName  = NULL;
  char *bedFName  = NULL;
  FILE *bedFp = NULL;

  char *dbUser = "ensro";
  char *dbPass = NULL;
  int   dbPort = 3306;

  char *dbHost = "ens-livemirror.internal.sanger.ac.uk";
  char *dbName = "homo_sapiens_core_71_37";

  char *assName = "GRCh37";

  char *chrName = "GL000206.1";


  int flags = 0;
  int   threads  = 1;

  initEnsC(argc, argv);

  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

// Ones without a val go here
    if (!strcmp(arg, "-U") || !strcmp(arg,"--ucsc_naming")) {
      flags |= M_UCSC_NAMING;
    } else {
// Ones with a val go in this block
      if (argNum == argc-1) {
        Bamcount_usage();
      }

      val = argv[++argNum];

      if (!strcmp(arg, "-i") || !strcmp(arg,"--in_file")) {
        StrUtil_copyString(&inFName,val,0);
      } else if (!strcmp(arg, "-b") || !strcmp(arg,"--bed")) {
        StrUtil_copyString(&bedFName,val,0);
      } else if (!strcmp(arg, "-h") || !strcmp(arg,"--host")) {
        StrUtil_copyString(&dbHost,val,0);
      } else if (!strcmp(arg, "-p") || !strcmp(arg,"--password")) {
        StrUtil_copyString(&dbPass,val,0);
      } else if (!strcmp(arg, "-P") || !strcmp(arg,"--port")) {
        dbPort = atoi(val);
      } else if (!strcmp(arg, "-n") || !strcmp(arg,"--name")) {
        StrUtil_copyString(&dbName,val,0);
      } else if (!strcmp(arg, "-u") || !strcmp(arg,"--user")) {
        StrUtil_copyString(&dbUser,val,0);
      } else if (!strcmp(arg, "-t") || !strcmp(arg,"--threads")) {
        threads = atoi(val);
      } else if (!strcmp(arg, "-a") || !strcmp(arg,"--assembly")) {
        StrUtil_copyString(&assName,val,0);
      } else if (!strcmp(arg, "-v") || !strcmp(arg,"--verbosity")) {
        verbosity = atoi(val);
// Temporary
      } else if (!strcmp(arg, "-c") || !strcmp(arg,"--chromosome")) {
        StrUtil_copyString(&chrName,val,0);
      } else {
        fprintf(stderr,"Error in command line at %s\n\n",arg);
        Bamcount_usage();
      }
    }
    argNum++;
  }

  if (verbosity > 0) {
    printf("Program for read counting of BAM reads against annotation \n"
           "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Mar 2013.\n");
  }

  if (!inFName) {
    Bamcount_usage();
  }

  if (bedFName) {
    bedFp = fopen(bedFName, "w");
  }
  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);

  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  slices = SliceAdaptor_fetchAll(sa, "toplevel", NULL, 0);
  nSlices = Vector_getNumElement(slices);

  if (Vector_getNumElement(slices) == 0) {
    fprintf(stderr, "Error: No slices.\n");
    exit(1);
  }


long long totalUsableReads = 1;

  fprintf(stderr, "Have %lld total usable reads\n",totalUsableReads);

  htsFile *in = hts_open(inFName, "rb");
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return 1;
  }

  hts_set_threads(in, threads);
  hts_idx_t *idx;
  idx = sam_index_load(in, inFName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM index file is not available.\n");
    return 1;
  }

  int i;
  for (i=0; i<Vector_getNumElement(slices); i++) {
    Slice *slice = Vector_getElementAt(slices,i);

    if (verbosity > 0) fprintf(stderr, "Working on '%s'\n",Slice_getSeqRegionName(slice));

    if (verbosity > 0) fprintf(stderr, "Stage 1 - retrieving annotation from database\n");
    Vector *genes = getGenes(slice, flags);

    IDHash *geneResultsHash = makeGeneResultsHash(genes);

    if (verbosity > 0) fprintf(stderr, "Stage 2 - counting reads\n");
    countReads(inFName, slice, in, idx, flags, genes, geneResultsHash, totalUsableReads, bedFp);
  }


  hts_idx_destroy(idx);
  hts_close(in);

  if (verbosity > 0) fprintf(stderr, "Done\n");
  if (bedFp) fclose(bedFp);
  return 0;
}

/*
 Program usage message
*/
void Bamcount_usage() {
  printf("bamcount \n"
         "  -i --in_file     Input BAM file to map from (string)\n"
         "  -U --ucsc_naming Input BAM file has 'chr' prefix on ALL seq region names (flag)\n"
         "  -h --host        Database host name for db containing mapping (string)\n"
         "  -n --name        Database name for db containing mapping (string)\n"
         "  -u --user        Database user (string)\n"
         "  -p --password    Database password (string)\n"
         "  -P --port        Database port (int)\n"
         "  -a --assembly    Assembly name (string)\n"
         "  -v --verbosity   Verbosity level (int)\n"
         "\n"
         "Notes:\n"
         "  -U will cause 'chr' to be prepended to all source seq_region names, except the special case MT which is changed to chrM.\n"
         "  -v Default verbosity level is 1. You can make it quieter by setting this to 0, or noisier by setting it > 1.\n"
         );
  exit(1);
}

int countReads(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags, Vector *origGenesVec, IDHash *geneResultsHash, long long countUsableReads, FILE *bedFp) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];
  char region_name[512];
  GeneResults *gr;
  Vector *genes = Vector_copy(origGenesVec);

  if (flags & M_UCSC_NAMING) {
    sprintf(region,"chr%s", Slice_getSeqRegionName(slice));
  } else {
    sprintf(region,"%s", Slice_getSeqRegionName(slice));
  }
  bam_hdr_t *header = bam_hdr_init();
  header = bam_hdr_read(in->fp.bgzf);
  ref = bam_name2id(header, region);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s\n", region);
    exit(1);
  }
  sprintf(region,"%s:%ld-%ld", region_name,
      Slice_getSeqRegionStart(slice),
      Slice_getSeqRegionEnd(slice));
  if (hts_parse_reg(region, &begRange, &endRange) == NULL) {
    fprintf(stderr, "Could not parse %s\n", region);
    exit(2);
  }
  bam_hdr_destroy(header);

  hts_itr_t *iter = sam_itr_queryi(idx, ref, begRange, endRange);
  bam1_t *b = bam_init1();

  long counter = 0;
  long overlapping = 0;
  long bad = 0;
  int startIndex = 0;
  while (bam_itr_next(in, iter, b) >= 0) {
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
      bad++;
      continue;
    }

    Vector *features = NULL;
    int end;
    end = bam_endpos(b);

    // There is a special case for reads which have zero length and start at begRange (so end at begRange ie. before the first base we're interested in).
    // That is the reason for the || end == begRange test
    if (end == begRange) {
      continue;
    }
    counter++;

    if (!(counter%1000000)) {
      if (verbosity > 1) { printf("."); }
      fflush(stdout);
    }

    int j;
    int done = 0;
    int hadOverlap = 0;

    for (j=startIndex; j < Vector_getNumElement(genes) && !done; j++) {
      Gene *gene = Vector_getElementAt(genes,j);
      if (!gene) {
        continue;
      }
      // Remember: b->core.pos is zero based!
      if (b->core.pos < Gene_getEnd(gene) && end >= Gene_getStart(gene)) {
        int k;

        int doneGene = 0;

        gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));
        Vector *geneExons = gr->exons;
        int m;
        for (m=0; m<Vector_getNumElement(geneExons) && !doneGene; m++) {
          Exon *exon = Vector_getElementAt(geneExons, m);

          if (b->core.pos < Exon_getEnd(exon) && end >= Exon_getStart(exon)) {
            if (features == NULL) {
              features = Bam_cigarToUngapped(b);
            }
            int n=0;
            for (n=0; n<Vector_getNumElement(features); n++) {
              DNAAlignFeature *daf = Vector_getElementAt(features, n);
              if (DNAAlignFeature_getStart(daf) <= Exon_getEnd(exon) &&
                  DNAAlignFeature_getEnd(daf) >= Exon_getStart(exon)) {
                Exon_setScore(exon, Exon_getScore(exon) + 1);
              }
            }
          }
        }
      } else if (Gene_getStart(gene) > end) {
        done = 1;
      } else if (Gene_getEnd(gene) < b->core.pos+1) {
        if (verbosity > 1) {
          printf("Removing gene %s (index %d) with extent %ld to %ld\n",
              Gene_getStableId(gene),
              gr->index,
              Gene_getStart(gene),
              Gene_getEnd(gene));
        }
        Vector_setElementAt(genes,j,NULL);

        // Magic (very important for speed) - move startIndex to first non null gene
        int n;
        startIndex = 0;
        for (n=0;n<Vector_getNumElement(genes);n++) {
          void *v = Vector_getElementAt(genes,n);

          if (v != NULL) {
            break;
          }
          startIndex++;
        }
        if (verbosity > 1) {
          printf("startIndex now %d\n",startIndex);
        }
      }
    }
    if (features) {
      // Temporary - use free to free DNAAlignFeature
      Vector_setFreeFunc(features, free);
      Vector_free(features);
    }
  }
  if (verbosity > 1) { printf("\n"); }

  // Print out read counts for what ever's left in the genes array
  int i;
  for (i=0;i<Vector_getNumElement(origGenesVec);i++) {
    Gene *gene = Vector_getElementAt(origGenesVec, i);

    gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));

    // Get max exon score
    double maxScore = 0.0;
    int j;
    for (j=0; j<Vector_getNumElement(gr->exons); j++) {
      Exon *exon = Vector_getElementAt(gr->exons, j);
      double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
      if (normExonScore > maxScore) {
        maxScore = normExonScore;
      }
    }
    fprintf(stderr, "Gene %s biotype %s num trans %d max normalised exon score %f\n", Gene_getStableId(gene),
        Gene_getBiotype(gene), Gene_getTranscriptCount(gene), maxScore);

    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *trans = Gene_getTranscriptAt(gene, j);
      Transcript_sort(trans);
      fprintf(stderr, "Transcript %s biotype %s num exon %d\n", Transcript_getStableId(trans),
          Transcript_getBiotype(trans), Transcript_getExonCount(trans));

      int k;
      double maxScoreInTrans = 0.0;
      for (k=0; k<Transcript_getExonCount(trans); k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
        if (normExonScore > maxScoreInTrans) {
          maxScoreInTrans = normExonScore;
        }
      }

      Exon *prevExon = NULL;
      Exon *nextExon = NULL;
      for (k=0; k<Transcript_getExonCount(trans); k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        if (k<Transcript_getExonCount(trans)-1) {
          nextExon = Transcript_getExonAt(trans, k+1);
        } else {
          nextExon = NULL;
        }

        double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
        fprintf(stderr, "Exon score for %s (%s) %s gsource %s %s %s tsource %s %s %s %ld %ld %d length %ld score %f score per base %f norm score/base %f",
            Gene_getStableId(gene),
            Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "",
            Gene_getBiotype(gene), Gene_getSource(gene),
            Transcript_getStableId(trans), Transcript_getBiotype(trans), Analysis_getLogicName(Transcript_getAnalysis(trans)),
            Exon_getStableId(exon), Slice_getSeqRegionName(Exon_getSlice(exon)), Exon_getSeqRegionStart(exon), Exon_getSeqRegionEnd(exon), Exon_getSeqRegionStrand((SeqFeature*)exon),
            Exon_getLength(exon), Exon_getScore(exon), Exon_getScore(exon)/Exon_getLength(exon),
            maxScore != 0.0 ? (((normExonScore)/maxScore)*1000.0) : 0.0);

        if (normExonScore < maxScoreInTrans/15.0 && prevExon!=NULL && nextExon!=NULL && Exon_getLength(exon) > 50) {
          double prevNormExonScore = Exon_getScore(prevExon)/Exon_getLength(prevExon);
          double nextNormExonScore = Exon_getScore(nextExon)/Exon_getLength(nextExon);
          if (prevNormExonScore/5.0 > normExonScore && nextNormExonScore/5.0 > normExonScore) {
            fprintf(stderr," LOWSCORE");
          }
        }
        fprintf(stderr,"\n");

        prevExon = exon;
      }
    }

    sam_itr_destroy(iter);
    bam_destroy1(b);
  }
}

Vector *Bam_cigarToUngapped(bam1_t *b) {
  int cigInd;
  int refPos;
  int readPos;
  uint32_t *cigar = bam_get_cigar(b);

 Vector *features = Vector_new();

  for (cigInd = readPos = 0, refPos = b->core.pos; cigInd < b->core.n_cigar; ++cigInd) {
    int k;
    int lenCigBlock = cigar[cigInd]>>4;
    int op          = cigar[cigInd]&0xf;

    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      DNAAlignFeature *daf = DNAAlignFeature_new();

      DNAAlignFeature_setStart(daf, refPos+1);
      DNAAlignFeature_setEnd(daf, refPos+lenCigBlock+1);

      Vector_addElement(features, daf);

      refPos += lenCigBlock; readPos += lenCigBlock;
    } else if (op == BAM_CDEL) {
      refPos += lenCigBlock;
    } else if (op == BAM_CSOFT_CLIP) {
      readPos += lenCigBlock;
    } else if (op == BAM_CHARD_CLIP) {
    } else if (op == BAM_CINS) {
       readPos += lenCigBlock;
    } else if (op == BAM_CREF_SKIP) {
       refPos += lenCigBlock;
    }
  }

  return features;
}
