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

#include "htslib/sam.h"
#define MY_FUSEDFLAG 32768
/*
 print out a bam1_t entry, particularly the flags (for debugging)
*/
void printBam(FILE *fp, bam1_t *b, bam_hdr_t *header) {
  fprintf(fp, "%s %s %d %d %s (%d) %d %d %d\t\tP %d PP %d U %d MU %d R %d MR %d R1 %d R2 %d S %d QC %d D %d U %d\n",
                                  bam_get_qname(b),
                                  header->target_name[b->core.tid],
                                  b->core.pos,
                                  bam_endpos(b),
                                  header->target_name[b->core.mtid],
                                  b->core.mtid,
                                  b->core.mpos,
                                  b->core.isize,
                                  bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)),
                                  b->core.flag & BAM_FPAIRED,
                                  b->core.flag & BAM_FPROPER_PAIR ? 1 : 0,
                                  b->core.flag & BAM_FUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FMUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FMREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FREAD1 ? 1 : 0,
                                  b->core.flag & BAM_FREAD2 ? 1 : 0,
                                  b->core.flag & BAM_FSECONDARY ? 1 : 0,
                                  b->core.flag & BAM_FQCFAIL ? 1 : 0,
                                  b->core.flag & BAM_FDUP ? 1 : 0,
                                  b->core.flag & MY_FUSEDFLAG ? 1 : 0
                                  );
  fflush(fp);
}

typedef struct GeneResultsStruct {
  int   index;
  long  score;
  Gene *gene;
  Vector *flatFeatures;
  Vector *exons;
  long  flatLength;
} GeneResults;

long long countReadsInFile(char *inFName) {
  htsFile *in = hts_open(inFName, "rb");
  long long nUsableReads = 0;

  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    exit(1);
  }

  bam1_t *b = bam_init1();

  int cnt = 0;
  while (bam_read1(in->fp.bgzf, b) > 0 && b->core.tid >= 0) {
    if (!(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) {
      nUsableReads++;
    }
    cnt++;
    if (!(cnt%1000000)) {
      fprintf(stderr,".");
      fflush(stdout);
    }
  }
  fprintf(stderr, "\n");

  hts_close(in);

  return nUsableReads;
}

int geneStartCompFunc(const void *one, const void *two) {
  Gene *g1 = *((Gene**)one);
  Gene *g2 = *((Gene**)two);

  return Gene_getStart(g1) - Gene_getStart(g2);
}

Vector *getGenes(Slice *slice, int flags) {
  Vector *genes;
  genes = Slice_getAllGenes(slice, NULL, NULL, 1, NULL, NULL);

  Vector_sort(genes, geneStartCompFunc);

  return genes;
}

Vector *flattenGene(Gene *gene) {
  SeqFeature *curBlock = NULL;
  Vector *exons = Gene_getAllExons(gene);
  Vector *blockFeatures = Vector_new();

  Vector_sort(exons, SeqFeature_startCompFunc);

  int i;
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    long exStart = Exon_getStart(exon) <= 0 ? 1 : Exon_getStart(exon);
    long exEnd   = Exon_getEnd(exon) <= 0 ?   1 : Exon_getEnd(exon);

    if (curBlock && SeqFeature_getEnd(curBlock) >= exStart) {
      if (exEnd > SeqFeature_getEnd(curBlock)) SeqFeature_setEnd(curBlock, exEnd);
    } else {
      curBlock = SeqFeature_new();
      SeqFeature_setStart(curBlock, exStart);
      SeqFeature_setEnd(curBlock, exEnd);
      Vector_addElement(blockFeatures, curBlock);
    }
  }

  Vector_sort(blockFeatures, SeqFeature_startCompFunc);
  Vector_free(exons);


  return blockFeatures;
}

// Pregenerate the hash for gene scores, so we don't have to test for existence within the BAM read loop
// Also means all genes will have a score hash entry
IDHash *makeGeneResultsHash(Vector *genes) {
  IDHash *geneResultsHash = IDHash_new(IDHASH_LARGE);

  int i;
  for (i=0; i < Vector_getNumElement(genes); i++) {
    Gene *gene = Vector_getElementAt(genes,i);
    GeneResults *gr;

    if ((gr = (GeneResults *)calloc(1,sizeof(GeneResults))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating GeneResults\n");
      exit(1);
    }
    gr->index = i;
    gr->gene  = gene;
    gr->flatFeatures = flattenGene(gene);
    int j=0;
    gr->flatLength = 0;
    for (j=0; j<Vector_getNumElement(gr->flatFeatures); j++) {
      SeqFeature *sf = Vector_getElementAt(gr->flatFeatures,j);
      gr->flatLength += SeqFeature_getLength(sf);
    }

    IDHash_add(geneResultsHash,Gene_getDbID(gene),gr);
  }

  return geneResultsHash;
}
