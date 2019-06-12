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

#include "RoughTranscripts.h"
#include <stdio.h>

#include <sysexits.h>
#include "ensc/EnsC.h"

#include "ensc/DBAdaptor.h"
#include "ensc/AnalysisAdaptor.h"
#include "ensc/GeneAdaptor.h"
#include "ensc/Vector.h"
#include "ensc/Slice.h"
#include "ensc/SliceAdaptor.h"
#include "ensc/StrUtil.h"
#include "ensc/DNAAlignFeature.h"
#include "ensc/Gene.h"
#include "ensc/Transcript.h"
#include "ensc/Attribute.h"
#include "ensc/MetaContainer.h"

#include "libconfig.h"
#ifdef HAVE_LIBTCMALLOC
#include "gperftools/tcmalloc.h"
#endif

#include "htslib/sam.h"
#include "htslib/hts.h"
/*
=head1 DESCRIPTION

This module takes intron spanning dna_align_features and combines them with
rough transcript models to build refined genes with CDS. The module produces
transcripts representing all possible combinations of introns and exons which
are then filtered according to the specifications in the config.
The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::RoughTranscripts
*/

#define EXIT_MEMORY 11
#define NO_GENES_WRITTEN 42

#define RSGVERSION "0.3.6-91"
#define RSGGENE_KEEP 16
#define RSG_DUPLICATE 32
#define RSG_PROCESSED 64
int dumpGenes(Vector *genes, int withSupport);

void RoughTranscripts_usage(int exit_code) {
  printf("RoughTranscripts \n"
         "  -c --config_file RoughTranscripts configuration file to read from\n"
         "  -i --input_id    Input id (slice name) to run on eg. chromosome:Oar_v3.1:17\n"
         "  -l --logic_name  Logic name for analysis block to run from configuration file\n"
         "  -d --dry_run     If specified, don't write to output db\n"
         "  -u --ucsc_naming If specified, add chr the name of sequence\n"
         "  -t --threads     Number of threads to use when reading the BAM files, default is 1\n"
         "  -v --verbosity   Verbosity level (int)\n"
         "  -V --version     Version of RoughTranscripts\n"
         "  -h --help        This help...\n"
         "\n"
         "Notes:\n"
         "  If no genes were written in the database, exit code is %i\n",
         NO_GENES_WRITTEN
         );
  exit(exit_code);
}

void RoughTranscripts_version() {
  printf("Version: %s\n", RSGVERSION);
  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  if (argc == 1) RoughTranscripts_usage(EX_USAGE);
  initEnsC(argc, argv);

  char *logicName;
  char *configFile;
  char *inputId;
  int   dryRun     = 0;
  int   threads  = 1;
  int   verbosity  = 0;
  int   ucsc_naming = 0;
  float min_span = 1.5;
  int   min_single_exon_length = 1000;
  int   min_length = 300;
  char *bam_file;
  int   dbport = 3306;
  char *dbhost;
  char *dbuser;
  char *dbpass;
  char *dbname;


  int argNum = 1;
  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

    // Ones without a val go here
    if (!strcmp(arg, "-V") || !strcmp(arg,"--version")) {
      RoughTranscripts_version();
    }
    else if (!strcmp(arg, "-h") || !strcmp(arg,"--help")) {
      RoughTranscripts_usage(EXIT_SUCCESS);
    }
    else if (!strcmp(arg, "-d") || !strcmp(arg,"--dry_run")) {
      dryRun = 1;
    }
    else if (!strcmp(arg, "-u") || !strcmp(arg,"--ucsc_naming")) {
      ucsc_naming = 1;
    }
    else {
      // Ones with a val go in this block
      if (argNum == argc-1) {
        fprintf(stderr, "Error: Expected a value after last command line argument\n");
        RoughTranscripts_usage(EX_USAGE);
      }

      val = argv[++argNum];

      if (!strcmp(arg, "-c") || !strcmp(arg,"--config_file")) {
        StrUtil_copyString(&configFile,val,0);
      }
      else if (!strcmp(arg, "-i") || !strcmp(arg,"--input_id")) {
        StrUtil_copyString(&inputId,val,0);
      }
      else if (!strcmp(arg, "-l") || !strcmp(arg,"--logic_name")) {
        StrUtil_copyString(&logicName,val,0);
      }
      else if (!strcmp(arg, "-t") || !strcmp(arg,"--threads")) {
        threads = atoi(val);
      }
      else if (!strcmp(arg, "-v") || !strcmp(arg,"--verbosity")) {
        verbosity = atoi(val);
      }
      else if (!strcmp(arg,"--dbhost")) {
        StrUtil_copyString(&dbhost, val, 0);
      }
      else if (!strcmp(arg,"--dbport")) {
        dbport = atoi(val);
      }
      else if (!strcmp(arg,"--dbuser")) {
        StrUtil_copyString(&dbuser,val,0);
      }
      else if (!strcmp(arg,"--dbpass")) {
        StrUtil_copyString(&dbpass,val,0);
      }
      else if (!strcmp(arg,"--dbname")) {
        StrUtil_copyString(&dbname,val,0);
      }
      else if (!strcmp(arg,"--bamfile")) {
        StrUtil_copyString(&bam_file,val,0);
      }
      else if (!strcmp(arg,"--min_single_exon_length")) {
        min_single_exon_length = atoi(val);
      }
      else if (!strcmp(arg,"--min_length")) {
        min_length = atoi(val);
      }
      else if (!strcmp(arg,"--min_span")) {
        min_span = atof(val);
      }
      else {
        fprintf(stderr,"Error in command line at %s\n\n",arg);
        RoughTranscripts_usage(EX_USAGE);
      }
    }

    argNum++;
  }
  int exit_code = EXIT_SUCCESS;

  RoughTranscripts *rsg = RoughTranscripts_new(configFile, logicName, inputId, threads, ucsc_naming, verbosity);
  RoughTranscripts_fetchInput(rsg);
  RoughTranscripts_run(rsg);

  if (RoughTranscripts_getOutput(rsg)) {
    if (verbosity > 0) {
      dumpGenes(RoughTranscripts_getOutput(rsg), 1);
    }
    if (dryRun) {
      fprintf(stderr,"DRY RUN mode - NOT writing genes to output db\n");
    }
    else {
      RoughTranscripts_writeOutput(rsg);
    }
#ifdef HAVE_LIBTCMALLOC
    if (verbosity > 0) tc_malloc_stats();
#endif
    Vector_setFreeFunc(rsg->output, Gene_free);
    Vector_free(rsg->output);
  }
  else {
    exit_code = NO_GENES_WRITTEN;
  }

#ifdef HAVE_LIBTCMALLOC
  if (verbosity > 0) tc_malloc_stats();
#endif
  return exit_code;
}

void RoughTranscripts_dumpOutput(RoughTranscripts *rsg) {
  Vector *output = RoughTranscripts_getOutput(rsg);

  Vector_sort(output, SeqFeature_startCompFunc);

  int i;
  for (i=0;i<Vector_getNumElement(output);i++) {
    Gene *gene = Vector_getElementAt(output, i);
    fprintf(stderr, "Gene\t%ld\t%ld\t%d\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
    int j;
    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *trans = Gene_getTranscriptAt(gene, j);
      fprintf(stderr, "Transcript\t%ld\t%ld\t%d\n", Transcript_getStart(trans), Transcript_getEnd(trans), Transcript_getStrand(trans));
      int k;
      for (k=0;k<Transcript_getExonCount(trans);k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        fprintf(stderr, "Exon\t%ld\t%ld\t%d\n", Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon));
      }
    }
    fprintf(stderr,"\n");
  }
}

void RoughTranscripts_initSetFuncs(RoughTranscripts *rsg) {
  rsg->funcHash = StringHash_new(STRINGHASH_SMALL);

  StringHash_add(rsg->funcHash, "BAM_FILE", SetFuncData_new(RoughTranscripts_setBamFile, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "TARGET_DB", SetFuncData_new(RoughTranscripts_setTargetDb, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "MIN_SPAN", SetFuncData_new(RoughTranscripts_setMinSpan, CONFIG_TYPE_FLOAT));
  StringHash_add(rsg->funcHash, "MAX_INTRON_SIZE", SetFuncData_new(RoughTranscripts_setMaxIntronSize, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MIN_SINGLE_EXON", SetFuncData_new(RoughTranscripts_setMinSingleExonLength, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "PAIRED", SetFuncData_new(RoughTranscripts_setDataPaired, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "STRANDED", SetFuncData_new(RoughTranscripts_setDataStranded, CONFIG_TYPE_INT));
}

void RoughTranscripts_dumpConfig(RoughTranscripts *rsg) {
  fprintf(stderr, "TARGET_DB\t\t%s\n", RoughTranscripts_getTargetDb(rsg));
  fprintf(stderr, "BAM_FILE\t\t%s\n", RoughTranscripts_getBamFile(rsg));
  fprintf(stderr, "MIN_SPAN\t\t%d\n", RoughTranscripts_getMinSpan(rsg));
  fprintf(stderr, "MAX_INTRON_SIZE\t\t%d\n", RoughTranscripts_getMaxIntronSize(rsg));
  fprintf(stderr, "MIN_SINGLE_EXON\t\t%d\n", RoughTranscripts_getMinSingleExonLength(rsg));
  fprintf(stderr, "PAIRED\t\t%d\n", RoughTranscripts_isDataPaired(rsg));
  fprintf(stderr, "STRANDED\t\t%d\n", RoughTranscripts_isDataStranded(rsg));
}

RoughTranscripts *RoughTranscripts_new(char *configFile, char *logicName, char *inputId, int threads, int ucsc_naming, int verbosity) {
  RoughTranscripts *rsg;

  if ((rsg = calloc(1,sizeof(RoughTranscripts))) == NULL) {
    fprintf(stderr,"Failed allocating RoughTranscripts\n");
    exit(EXIT_MEMORY);
  }

  RoughTranscripts_initSetFuncs(rsg);

  rsg->adaptorAliasHash = StringHash_new(STRINGHASH_SMALL);

  if (configFile) {
    RunnableDB_readAndCheckConfig(rsg, configFile, "Config.ROUGHTRANSCRIPTS_CONFIG_BY_LOGIC", logicName);
  } else {
    fprintf(stderr, "WARNING: Running without reading config (config reading not implemented)\n");
    exit(EX_DATAERR);
  }
  RoughTranscripts_setInputId(rsg, inputId);
  RoughTranscripts_setAnalysis(rsg, RoughTranscripts_createAnalysisObject(rsg, logicName));

  RoughTranscripts_setUcscNaming(rsg, ucsc_naming);
  RoughTranscripts_setThreads(rsg, threads);
  RoughTranscripts_setVerbosity(rsg, verbosity);

  return rsg;
}

void RunnableDB_readAndCheckConfig(RoughTranscripts *rsg, char *configFile, char *blockName, char *logicName) {
  config_t cfg;
  config_setting_t *cfgBlock;

  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if (!config_read_file(&cfg, configFile)) {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EX_CONFIG);
  }

  // For now specify location of Databases config in config file and read it here
  config_setting_t *databasesFileSetting = config_lookup(&cfg, "Config.DATABASES_FILE");
  if (databasesFileSetting == NULL) {
    fprintf(stderr,"Missing config setting DATABASES_FILE\n");
    exit(EX_CONFIG);
  }
  const char *databasesFile = config_setting_get_string(databasesFileSetting);
  RunnableDB_readDatabaseConfig(rsg, (char *)databasesFile);

  cfgBlock = config_lookup(&cfg, blockName);
  if (cfgBlock == NULL) {
    fprintf(stderr,"Missing config block %s\n", blockName);
    exit(EX_CONFIG);
  }


// HACK: For now use passed in logicName rather than doing through analysis
  Utilities_parseConfig(rsg, cfgBlock, logicName, 0 /*ignoreThrow*/);
}

/*
=head2 fetch_sequence

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string, name
  Arg [3]   : Bio::EnsEMBL::DBAdaptor
  Arg [4]   : arrayref of logic_name if sequence is to be masked
  Arg [5]   : Boolean for softmasking if sequence is to be softmasked
  Function  : gets sequence from specifed database
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: none
  Example   :

=cut
*/
// This should go in RunnableDB
Slice *RoughTranscripts_fetchSequence(RoughTranscripts *rsg, char *name, DBAdaptor *db, Vector *repeatMaskTypes, int softMask) {
  if (!db) {
    db = RoughTranscripts_getDb(rsg);
  }

  if (!name) {
    name = RoughTranscripts_getInputId(rsg);
  }


  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(db);
  Slice *slice = SliceAdaptor_fetchByName(sa, name);

  if (slice == NULL) {
    fprintf(stderr,"Failed to fetch slice %s\n",name);
    exit(EX_IOERR);
  }

  fprintf(stderr, "slice name = %s\n", Slice_getName(slice));

  if (repeatMaskTypes != NULL && Vector_getNumElement(repeatMaskTypes)) {
    fprintf(stderr,"RepeatMasked sequence fetching not implemented yet - bye\n");
    exit(EX_SOFTWARE);
  }

  return slice;
}

void RoughTranscripts_fetchInput(RoughTranscripts *rsg) {
  int verbosity = RoughTranscripts_getVerbosity(rsg);

  if (verbosity > 0) fprintf(stderr, "INFO: FETCH INPUT\n");
  RoughTranscripts_setDb(rsg, BaseGeneBuild_getDbAdaptor(rsg, "REFERENCE_DB", 0, 0));

  DBAdaptor *db = RoughTranscripts_getDb(rsg);

  Slice *slice = RoughTranscripts_fetchSequence(rsg, RoughTranscripts_getInputId(rsg), NULL, NULL, 0);
  RoughTranscripts_setSlice(rsg, slice);

  RoughTranscripts_bamToExonFeatures(rsg);
}


void RoughTranscripts_run (RoughTranscripts *rsg) {
  int verbosity = RoughTranscripts_getVerbosity(rsg);
  int i;
  int j;
  if (verbosity > 0) fprintf(stderr, "INFO: RUN\n");
  Vector *forward = RoughTranscripts_getForwardExons(rsg);
  Vector *reverse = RoughTranscripts_getReverseExons(rsg);
  if (Vector_getNumElement(forward)) {
    if (verbosity > 0) fprintf(stderr, "INFO: PROCESS FORWARD %d\n", Vector_getNumElement(forward));
    RoughTranscripts_processClusters(rsg, forward);
    Vector_free(forward);
  }
  if (Vector_getNumElement(reverse)) {
    if (verbosity > 0) fprintf(stderr, "INFO: PROCESS REVERSE %d\n", Vector_getNumElement(reverse));
    for (i = 0; i < Vector_getNumElement(reverse); i++) {
      fprintf(stderr, "%5d %10ld %10ld\n", i, DAFCluster_getStart(Vector_getElementAt(reverse, i)), DAFCluster_getEnd(Vector_getElementAt(reverse, i)));
      DAFCluster *cluster = Vector_getElementAt(reverse, i);
      for (j = 0; j < Vector_getNumElement(cluster->links); j++) {
        fprintf(stderr, "%7d %10ld %10ld\n", j, DAFCluster_getStart(Vector_getElementAt(cluster->links, j)), DAFCluster_getEnd(Vector_getElementAt(cluster->links, j)));
      }
    }
    RoughTranscripts_processClusters(rsg, reverse);
    Vector_free(reverse);
  }
}


void RoughTranscripts_writeOutput(RoughTranscripts *rsg) {
  int verbosity = RoughTranscripts_getVerbosity(rsg);
  if (verbosity > 0) fprintf(stderr, "INFO: WRITE OUTPUT\n");
  Vector *output = RoughTranscripts_getOutput(rsg);
  if (Vector_getNumElement(output)) {
    DBAdaptor *outdb = BaseGeneBuild_getDbAdaptor(rsg, RoughTranscripts_getTargetDb(rsg), 0, 0);
    GeneAdaptor *geneAdaptor = DBAdaptor_getGeneAdaptor(outdb);

    Analysis *anal = RoughTranscripts_getAnalysis(rsg);
    if (verbosity > 0) fprintf(stderr,"Analysis set to %s\n", Analysis_getLogicName(anal));

    int i;
    for (i=0; i<Vector_getNumElement(output); i++) {
      Gene *gene = Vector_getElementAt(output, i);
      GeneAdaptor_store(geneAdaptor, gene, 1);
    }
  }
  else {
    if (verbosity) fprintf(stderr, "No genes to write");
    exit(NO_GENES_WRITTEN);
  }
}

void RoughTranscripts_addToOutput(RoughTranscripts *rsg, Gene *gene) {
  if (rsg->output == NULL) {
    rsg->output = Vector_new();
  }
  gene->flags |= RSGGENE_KEEP;
  Vector_addElement(rsg->output, gene);
}

Vector *RoughTranscripts_getOutput(RoughTranscripts *rsg) {
  return rsg->output;
}

void RoughTranscripts_setAnalysis(RoughTranscripts *rsg, Analysis *analysis) {
  rsg->analysis = analysis;
}

Analysis *RoughTranscripts_getAnalysis(RoughTranscripts *rsg) {
  return rsg->analysis;
}

Analysis *RoughTranscripts_createAnalysisObject(RoughTranscripts *rsg, char *logicName) {
  DBAdaptor *outDb = BaseGeneBuild_getDbAdaptor(rsg, RoughTranscripts_getTargetDb(rsg), 0, 0);

  AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(outDb);
  Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa, logicName);

  if (analysis != NULL) {
     return analysis;
  }
  // need to create analysis object first
  analysis = Analysis_new();
  char *chP;
  Analysis_setLogicName(analysis, StrUtil_copyString(&chP, logicName, 0));
  Analysis_setModule(analysis, StrUtil_copyString(&chP, "RoughTranscripts", 0));

  return analysis;
}


void RoughTranscripts_bamToExonFeatures(RoughTranscripts *rsg) {
  int verbosity = RoughTranscripts_getVerbosity(rsg);
  Slice *slice = RoughTranscripts_getSlice(rsg);
  Vector *current;
  Vector *forward = Vector_new();
  Vector *reverse = Vector_new();
  int cluster_count = 1;
  int read_count = 0;
  int total_read_count = 0;
  int realStrand = -1;
  int paired = RoughTranscripts_isDataPaired(rsg);
  int stranded_reads = RoughTranscripts_isDataStranded(rsg);
  Analysis *analysis = RoughTranscripts_getAnalysis(rsg);
  int lastPairing = 0;
  int *pLastPairing = &lastPairing;
  int first_mate_strand = -1;

  char *intronFile = RoughTranscripts_getBamFile(rsg);
  char region[2048];
  char region_name[1024];
  int ref;
  int begRange;
  int endRange;
  int i = 0;
  int start = 0;
  int end = 0;
  int createNewCluster = 1;

  htsFile *sam = hts_open(intronFile, "rb");
  if (sam == NULL) {
    fprintf(stderr, "Bam file %s not found\n", intronFile);
    exit(EX_NOINPUT);
  }
  if (verbosity > 0) fprintf(stderr,"Opened bam file %s\n", intronFile);

  hts_set_threads(sam, RoughTranscripts_getThreads(rsg));
  if (verbosity > 0) fprintf(stderr,"Setting number of threads to %d\n", RoughTranscripts_getThreads(rsg));
  hts_idx_t *idx;
  idx = sam_index_load(sam, intronFile); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM index file is not available.\n");
    exit(EX_NOINPUT);
  }
  if (verbosity > 0) fprintf(stderr,"Opened bam index for %s\n", intronFile);

  if (RoughTranscripts_getUcscNaming(rsg) == 0) {
    sprintf(region_name,"%s", Slice_getSeqRegionName(slice));
  } else {
    sprintf(region_name,"chr%s", Slice_getSeqRegionName(slice));
  }
  bam_hdr_t *header = bam_hdr_init();
  header = bam_hdr_read(sam->fp.bgzf);
  ref = bam_name2id(header, region_name);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s\n", region_name);
    exit(EX_DATAERR);
  }
  sprintf(region,"%s:%ld-%ld", region_name, Slice_getSeqRegionStart(slice),
      Slice_getSeqRegionEnd(slice));

  if (hts_parse_reg(region, &begRange, &endRange) == NULL) {
    fprintf(stderr, "Could not parse %s\n", region);
    exit(EX_DATAERR);
  }
  if (verbosity > 0) fprintf(stderr,"Parsed region for region %s\n", region);
  hts_itr_t *iter = sam_itr_queryi(idx, ref, begRange, endRange);
  bam1_t *read = bam_init1();

  if (verbosity > 0) fprintf(stderr,"before bam read loop\n");
  while (bam_itr_next(sam, iter, read) >= 0) {
    ++total_read_count;
    // Skip unmapped read
    if ((read->core.flag & BAM_FUNMAP) != 0) {
      continue;
    }
    // ignore unspliced reads if the bam file is a mixture of spliced and
    // unspliced reads
    if (bam_aux_get(read, "XS")) {
      continue;
    }
    ++read_count;

    if (stranded_reads) {
      realStrand = bam_is_rev(read) == 1 ? -1 : 1;
      if ((read->core.flag & BAM_FREAD1) != 0) {
        realStrand *= first_mate_strand;
      }
    }
    start = read->core.pos+1;
    end = bam_endpos(read)+1;
    current = realStrand == 1 ? forward : reverse;
    createNewCluster = 1;
    for (i = Vector_getNumElement(current)-1; i > -1; i--) {
      DAFCluster *cluster = Vector_getElementAt(current, i);
      if (DAFCluster_getStart(cluster) <= end && DAFCluster_getEnd(cluster) >= start) {
        if (DAFCluster_getEnd(cluster) < end) {
          DAFCluster_setEnd(cluster, end);
        }
        DAFCluster_incrementScore(cluster);
        if (paired && read->core.isize < 0) {
          RoughTranscripts_processMatePairs(read, cluster, current, i, pLastPairing);
        }
        createNewCluster = 0;
        break;
      }
      else if (DAFCluster_getEnd(cluster) < start) {
        break;
      }
    }
    if (createNewCluster) {
      char hitSeqName[32];
      sprintf(hitSeqName, "%i", cluster_count++);
      DNAAlignFeature *feat = DNAAlignFeature_new();
      DNAAlignFeature_setSlice     (feat, slice);
      DNAAlignFeature_setStart     (feat, start);
      DNAAlignFeature_setEnd       (feat, end);
      DNAAlignFeature_setStrand    (feat, realStrand);
      DNAAlignFeature_setHitSeqName(feat, hitSeqName);
      DNAAlignFeature_setHitStart  (feat, 1);
      DNAAlignFeature_setHitStrand (feat, 1);
      DNAAlignFeature_setHitEnd    (feat, read->core.l_qseq);
      DNAAlignFeature_setAnalysis  (feat, analysis);
      DNAAlignFeature_setScore     (feat, 1);
      DNAAlignFeature_sethCoverage (feat, 0);
      DAFCluster *cluster = DAFCluster_new(feat);
      Vector_addElement(current, cluster);
      if (paired && read->core.isize < 0) {
        RoughTranscripts_processMatePairs(read, cluster, current, Vector_getNumElement(current)-1, pLastPairing);
      }
    }
  }
  if (verbosity > 0) fprintf(stderr,"after bam read loop\n");
  bam_destroy1(read);
  sam_itr_destroy(iter);

  hts_idx_destroy(idx);
  bam_hdr_destroy(header);
  hts_close(sam);

  RoughTranscripts_setForwardExons(rsg, forward);
  RoughTranscripts_setReverseExons(rsg, reverse);

  if (verbosity > 0) {
    fprintf(stderr, "Processed %d reads out of %d total reads\n", read_count, total_read_count);
    fprintf(stderr, "Got %d exons on the forward strand\n", Vector_getNumElement(forward));
    fprintf(stderr, "Got %d exons on the reverse strand\n", Vector_getNumElement(reverse));
  }
}


void RoughTranscripts_processMatePairs(bam1_t *read, DAFCluster *cluster, Vector *current, int index, int *last_pairing) {
  int mate_start = read->core.mpos+1;
  int mate_end = read->core.mpos+2;
  int j = 0;
  int k = 0;
  int low_boundary = 0;
  int high_boundary = 0;
  int start = 0;
  int end = 0;
  int isize = 0;
  DAFCluster *mate_cluster;
  if (!(mate_end >= DAFCluster_getStart(cluster) && mate_start <= DAFCluster_getEnd(cluster))) {
    j = index-1;
    if (DAFCluster_getStart(Vector_getElementAt(current, *last_pairing)) <= mate_end
        && DAFCluster_getEnd(Vector_getElementAt(current, *last_pairing)) >= mate_start) {
      mate_cluster = Vector_getElementAt(current, *last_pairing);
      if (Vector_getNumElement(mate_cluster->links)) {
        for (k = Vector_getNumElement(mate_cluster->links)-1; k > -1; k--) {
          if (Vector_getElementAt(mate_cluster->links, k) == cluster) {
            return;
          }
        }
      }
      Vector_addElement(mate_cluster->links, cluster);
      return;
    }
    else if (DAFCluster_getStart(Vector_getElementAt(current, j)) <= mate_end
        && DAFCluster_getEnd(Vector_getElementAt(current, j)) >= mate_start) {
      mate_cluster = Vector_getElementAt(current, j);
      *last_pairing = j;
      if (Vector_getNumElement(mate_cluster->links)) {
        for (k = Vector_getNumElement(mate_cluster->links)-1; k > -1; k--) {
          if (Vector_getElementAt(mate_cluster->links, k) == cluster) {
            return;
          }
        }
      }
      Vector_addElement(mate_cluster->links, cluster);
      return;
    }
    else {
      low_boundary = 0;
      high_boundary = j;
      start = read->core.pos+1;
      isize = read->core.isize;
      j = (mate_start/DAFCluster_getStart(Vector_getElementAt(current, j)))*high_boundary;
      while (!(DAFCluster_getStart(Vector_getElementAt(current, j)) <= mate_end && DAFCluster_getEnd(Vector_getElementAt(current, j)) >= mate_start)) {
        if (DAFCluster_getStart(Vector_getElementAt(current, j)) > mate_end) {
          high_boundary = j;
          j = high_boundary-((mate_start/DAFCluster_getStart(Vector_getElementAt(current, j)))*(high_boundary-low_boundary) || 1);
        }
        else {
          low_boundary = j;
          j = low_boundary+((mate_start/DAFCluster_getStart(Vector_getElementAt(current, high_boundary)))*(high_boundary-low_boundary) || 1);
        }
      }
      mate_cluster = Vector_getElementAt(current, j);
      *last_pairing = j;
      if (Vector_getNumElement(mate_cluster->links)) {
        for (k = Vector_getNumElement(mate_cluster->links)-1; k > -1; k--) {
          if (Vector_getElementAt(mate_cluster->links, k) == cluster) {
            return;
          }
        }
      }
      Vector_addElement(mate_cluster->links, cluster);
      return;
    }
  }
}

void RoughTranscripts_processClusters (RoughTranscripts *rsg, Vector *clusters) {
  Vector *transcripts = Vector_new();
  int cluster_index = Vector_getNumElement(transcripts);
  int i = 0;
  int j = 0;
  int last_end = 0;
  int tmp_index = -1;
  int max_intron_length = RoughTranscripts_getMaxIntronSize(rsg);
  int min_length = RoughTranscripts_getMinLength(rsg);
  int min_single_exon_length = RoughTranscripts_getMinSingleExonLength(rsg);
  float min_span = RoughTranscripts_getMinSpan(rsg);
  if (RoughTranscripts_isDataPaired(rsg)) {
    for (i = 0; i < Vector_getNumElement(clusters); i++) {
      DAFCluster *cluster = Vector_getElementAt(clusters, i);
      if (DNAAlignFeature_gethCoverage(cluster->feature) == 0) {
        if (Vector_getNumElement(cluster->links)) {
          cluster_index = Vector_getNumElement(transcripts);
          tmp_index = RoughTranscripts_processLink(cluster, transcripts, cluster_index);
          if (tmp_index != cluster_index) {
            exit(EX_DATAERR);
          }
        }
        else {
          DAFCluster_incrementHCoverage(cluster);
          Vector *transcript = Vector_new();
          Vector_addElement(transcript, cluster);
          Vector_addElement(transcripts, transcript);
          cluster_index = Vector_getNumElement(transcripts);
        }
      }
    }
  }
  else {
    last_end = DAFCluster_getEnd(Vector_getElementAt(clusters, 0));
    Vector_addElement(transcripts, Vector_new());
    Vector_addElement(Vector_getElementAt(transcripts, cluster_index), Vector_getElementAt(clusters, 0));
    for (i = 1; i < Vector_getNumElement(clusters); i++) {
      DAFCluster *cluster = Vector_getElementAt(clusters, i);
      if (DAFCluster_getStart(cluster)-last_end > max_intron_length) {
        Vector_addElement(transcripts, Vector_new());
        ++cluster_index;
      }
      last_end = DAFCluster_getEnd(cluster);
      Vector_addElement(Vector_getElementAt(transcripts, cluster_index), cluster);
    }
  }
//  for (i = 0; i < Vector_getNumElement(clusters); i++) {
//    DNAAlignFeature_sethCoverage(((DAFCluster *)Vector_getElementAt(clusters, i))->feature, 0);
//  }
  for (i = 0; i < Vector_getNumElement(transcripts); i++) {
    Vector *exons = Vector_getElementAt(transcripts, i);
    if (Vector_getNumElement(exons)) {
      Transcript *transcript = Transcript_new();
      for (j = 0; j < Vector_getNumElement(exons);j++) {
        DAFCluster *cluster = Vector_getElementAt(exons, j);
        DNAAlignFeature *feature = cluster->feature;
        Exon *exon = ExonUtils_createExon(
          DNAAlignFeature_getStart(feature),
          DNAAlignFeature_getEnd(feature),
          -1,
          -1,
          DNAAlignFeature_getStrand(feature),
          DNAAlignFeature_getAnalysis(feature),
          NULL,
          0,
          DNAAlignFeature_getSlice(feature),
          NULL,
          0
        );
        char cigar_string[32];
        sprintf(cigar_string, "%dM", Exon_getLength(exon));
        DNAAlignFeature_setCigarString(feature, cigar_string);
        Vector *tmpSf = Vector_new();
        Vector_addElement(tmpSf, feature);
        Exon_addSupportingFeatures(exon, tmpSf);
        Vector_free(tmpSf);
        Transcript_addExon(transcript, exon, 0);
        DAFCluster_free(cluster);
      }
      if (RoughTranscripts_isTranscriptValid(transcript, min_length, min_single_exon_length, min_span)) {
        TranscriptUtils_addExonPadding(transcript);
        Transcript_setAnalysis(transcript, Exon_getAnalysis((Exon *)Transcript_getExonAt(transcript, 0)));
        Transcript_setBiotype(transcript, Analysis_getLogicName(Transcript_getAnalysis(transcript)));
        Gene *gene = Gene_new();
        Gene_addTranscript(gene, transcript);
        Gene_setAnalysis(gene, Transcript_getAnalysis(transcript));
        Gene_setBiotype(gene, Transcript_getBiotype(transcript));
        Gene_setSource(gene, "ensembl");
        RoughTranscripts_addToOutput(rsg, gene);
      }
      else {
        fprintf(stderr, "WARNING: Transcript not valid Start: %ld End: %ld Strand: %d Num Exons: %d Num Clusters: %d\n", Transcript_getStart(transcript), Transcript_getEnd(transcript), Transcript_getStrand(transcript), Transcript_getExonCount(transcript), Vector_getNumElement(exons));
        Transcript_free(transcript);
      }
    }
    Vector_free(exons);
  }
  Vector_free(transcripts);
}

int RoughTranscripts_isTranscriptValid (Transcript *transcript, int min_length, int min_single_exon_length, float min_span) {
  if (Transcript_getLength(transcript) < min_length) {
    return 0;
  }
  if (Transcript_getExonCount(transcript) == 1) {
    if (Transcript_getLength(transcript) < min_single_exon_length) {
      return 0;
    }
  }
  else {
    // filter span on multiexon genes
    if ((float) ((Transcript_getEnd(transcript)-Transcript_getStart(transcript)+1)/Transcript_getLength(transcript)) < min_span) {
      if (Transcript_getLength(transcript) < min_single_exon_length) {
        return 0;
      }
    }
  }
  return 1;
}



int RoughTranscripts_processLink(DAFCluster *object, Vector *cluster, int index) {
  int i = 0;
  int j = 0;
  int found = 0;
  DAFCluster_incrementHCoverage(object);
  if (object->index != -1) {
    if (index != object->index) {
      if (Vector_getNumElement(cluster) == index) {
        Vector_addElement(cluster, Vector_new());
      }
      Vector *old_cluster = Vector_getElementAt(cluster, object->index);
      Vector *new_cluster = Vector_getElementAt(cluster, index);
      for (i = 0; i < Vector_getNumElement(old_cluster); i++) {
        DAFCluster *elm = Vector_getElementAt(old_cluster, i);
        found = 0;
        for (j = 0; j < Vector_getNumElement(new_cluster); j++) {
          if (Vector_getElementAt(new_cluster, j) == elm) {
            found = 1;
            break;
          }
        }
        if (! found) {
          Vector_addElement(new_cluster, Vector_getElementAt(old_cluster, i));
          elm->index = index;
        }
      }
      Vector_removeAll(old_cluster);
    }
  }
  else {
    object->index = index;
    for (i = 0; i < Vector_getNumElement(object->links); i++) {
      j = RoughTranscripts_processLink(Vector_getElementAt(object->links, i), cluster, index);
      if (j != index) {
        index = j;
      }
    }
    found = 0;
    if (Vector_getNumElement(cluster) != index) {
      Vector *new_cluster = Vector_getElementAt(cluster, index);
      for (i = 0; i < Vector_getNumElement(new_cluster); i++) {
        if (Vector_getElementAt(new_cluster, i) == object) {
          found = 1;
          break;
        }
      }
      if (! found) {
        Vector_addElement(new_cluster, object);
      }
    }
    else {
      Vector *new_cluster = Vector_new();
      Vector_addElement(new_cluster, object);
      Vector_addElement(cluster, new_cluster);
    }
  }
  return index;
}



//##################################################################
//# Containers



Slice *RoughTranscripts_getSlice(RoughTranscripts *rsg) {
  return rsg->chrSlice;
}

void RoughTranscripts_setSlice(RoughTranscripts *rsg, Slice *slice) {
  rsg->chrSlice = slice;
}


//####################################
//# config variable holders
//####################################

void RoughTranscripts_setTargetDb(RoughTranscripts *rsg, char *targetDb) {
  if (targetDb[0] != '\0') {
    if (rsg->targetDb != NULL) {
      free(rsg->targetDb);
    }
    rsg->targetDb = StrUtil_copyString(&rsg->targetDb, targetDb, 0);
  }
}

char *RoughTranscripts_getTargetDb(RoughTranscripts *rsg) {
  return rsg->targetDb;
}

void RoughTranscripts_setBamFile(RoughTranscripts *rsg, char *file) {
  if (file[0] != '\0') {
    if (rsg->bamFile != NULL) {
      free(rsg->bamFile);
    }
    rsg->bamFile = StrUtil_copyString(&rsg->bamFile, file, 0);
  }
}

char *RoughTranscripts_getBamFile(RoughTranscripts *rsg) {
  return rsg->bamFile;
}

void RoughTranscripts_setDataPaired(RoughTranscripts *rsg, int isPaired) {
  rsg->isPaired = isPaired;
}

int RoughTranscripts_isDataPaired(RoughTranscripts *rsg) {
  return rsg->isPaired;
}

void RoughTranscripts_setDataStranded(RoughTranscripts *rsg, int isStranded) {
  rsg->isStranded = isStranded;
}

int RoughTranscripts_isDataStranded(RoughTranscripts *rsg) {
  return rsg->isStranded;
}

void RoughTranscripts_setMinLength(RoughTranscripts *rsg, int minLength) {
  rsg->minLength = minLength;
}

int RoughTranscripts_getMinLength(RoughTranscripts *rsg) {
  return rsg->minLength;
}

void RoughTranscripts_setMinSpan(RoughTranscripts *rsg, float minSpan) {
  rsg->minSpan = minSpan;
}

float RoughTranscripts_getMinSpan(RoughTranscripts *rsg) {
  return rsg->minSpan;
}

void RoughTranscripts_setMaxIntronSize(RoughTranscripts *rsg, int maxIntronSize) {
  rsg->maxIntronSize = maxIntronSize;
}

int RoughTranscripts_getMaxIntronSize(RoughTranscripts *rsg) {
  return rsg->maxIntronSize;
}

void RoughTranscripts_setMinSingleExonLength(RoughTranscripts *rsg, int minSingleExonLength) {
  rsg->minSingleExonLength = minSingleExonLength;
}

int RoughTranscripts_getMinSingleExonLength(RoughTranscripts *rsg) {
  return rsg->minSingleExonLength;
}

void RoughTranscripts_setVerbosity(RoughTranscripts *rsg, int verbosity) {
  rsg->verbosity = verbosity;
}

int RoughTranscripts_getVerbosity(RoughTranscripts *rsg) {
  return rsg->verbosity;
}

void RoughTranscripts_setThreads(RoughTranscripts *rsg, int threads) {
  rsg->threads = threads;
}

int RoughTranscripts_getThreads(RoughTranscripts *rsg) {
  return rsg->threads;
}

void RoughTranscripts_setUcscNaming(RoughTranscripts *rsg, int ucsc_naming) {
  rsg->ucsc_naming = ucsc_naming;
}

int RoughTranscripts_getUcscNaming(RoughTranscripts *rsg) {
  return rsg->ucsc_naming;
}

void RoughTranscripts_setDb(RoughTranscripts *rsg, DBAdaptor *db) {
  rsg->db = db;
}

DBAdaptor *RoughTranscripts_getDb(RoughTranscripts *rsg) {
  return rsg->db;
}

void RoughTranscripts_setInputId(RoughTranscripts *rsg, char *inputId) {
    if (rsg->inputId != NULL) {
      free(rsg->inputId);
    }
  rsg->inputId = StrUtil_copyString(&rsg->inputId, inputId, 0);
}

char *RoughTranscripts_getInputId(RoughTranscripts *rsg) {
  return rsg->inputId;
}

void RoughTranscripts_setForwardExons(RoughTranscripts *rsg, Vector *exons) {
  rsg->forward = exons;
}

Vector *RoughTranscripts_getForwardExons(RoughTranscripts *rsg) {
  return rsg->forward;
}

void RoughTranscripts_setReverseExons(RoughTranscripts *rsg, Vector *exons) {
  rsg->reverse = exons;
}

Vector *RoughTranscripts_getReverseExons(RoughTranscripts *rsg) {
  return rsg->reverse;
}

Exon *ExonUtils_createExon(long start, long end, int phase, int endPhase, int strand, Analysis *analysis, Vector *supportingFeatures, IDType dbId, Slice *slice, char *stableId, int version) {

  Exon *newExon = Exon_new();

  Exon_setStart(newExon, start);
  Exon_setEnd(newExon, end);
  Exon_setPhase(newExon, phase);
  Exon_setEndPhase(newExon, endPhase);
  Exon_setStrand(newExon, strand);
  Exon_setAnalysis(newExon, analysis);
  Exon_setDbID(newExon, dbId);
  Exon_setSlice(newExon, slice);
  if (stableId) Exon_setStableId(newExon, stableId);
  if (version)  Exon_setVersion(newExon, version);

  if (supportingFeatures) {
    Exon_addSupportingFeatures(newExon, supportingFeatures);
  }

  return newExon;
}

void *TranscriptUtils_addExonPadding(Transcript *transcript) {
  int i = 0;
  int last_end = 0;
  int trim = 0;

  Vector *sortedExons = Vector_new();
  for (i = 0; i < Transcript_getExonCount(transcript); i++) {
    Vector_addElement(sortedExons, Transcript_getExonAt(transcript, i));
  }
  Vector_sort(sortedExons, SeqFeature_startCompFunc);
  for (i = 0; i < Vector_getNumElement(sortedExons); i++) {
    Exon_setStart((Exon *)Vector_getElementAt(sortedExons, i), Exon_getStart((Exon *)Vector_getElementAt(sortedExons, i))-20);
    Exon_setEnd((Exon *)Vector_getElementAt(sortedExons, i), Exon_getEnd((Exon *)Vector_getElementAt(sortedExons, i))+20);
    if (Exon_getStart((Exon *)Vector_getElementAt(sortedExons, i)) < 1) {
      Exon_setStart((Exon *)Vector_getElementAt(sortedExons, i), 1);
    }
    if (Exon_getEnd((Exon *)Vector_getElementAt(sortedExons, i)) > Slice_getSeqRegionLength(Transcript_getSlice(transcript))) {
      Exon_setEnd((Exon *)Vector_getElementAt(sortedExons, i), Slice_getSeqRegionLength(Transcript_getSlice(transcript)));
    }
    if (last_end) {
      if (last_end >= Exon_getStart((Exon *)Vector_getElementAt(sortedExons, i))) {
        trim = (last_end - Exon_getStart((Exon *)Vector_getElementAt(sortedExons, i)))/2;
        Exon_setStart((Exon *)Vector_getElementAt(sortedExons, i), Exon_getStart((Exon *)Vector_getElementAt(sortedExons, i))-trim+1);
        Exon_setEnd((Exon *)Vector_getElementAt(sortedExons, i-1), Exon_getEnd((Exon *)Vector_getElementAt(sortedExons, i-1))-1);
      }
    }
    last_end = Exon_getEnd((Exon *)Vector_getElementAt(sortedExons, i));
  }
  Vector_free(sortedExons);
}


int dumpGenes(Vector *genes, int withSupport) {
  FILE *fp = stderr;
  int i;
  int failed = 0;
  for (i=0;i<Vector_getNumElement(genes) && !failed;i++) {
    Gene *g = Vector_getElementAt(genes,i);
    fprintf(fp,"Gene %s (%s) type %s logicname %s coords: %ld %ld %d\n",Gene_getStableId(g),(Gene_getDisplayXref(g) ? DBEntry_getDisplayId(Gene_getDisplayXref(g)) : ""),Gene_getBiotype(g), Analysis_getLogicName(Gene_getAnalysis(g)), Gene_getStart(g),Gene_getEnd(g),Gene_getStrand(g));

    int j;
    for (j=0;j<Gene_getTranscriptCount(g);j++) {
      Transcript *t = Gene_getTranscriptAt(g,j);
      int k;

      fprintf(fp," Trans %s coords: %ld %ld %d biotype: %s logicname: %s\n",
              Transcript_getStableId(t),
              Transcript_getStart(t),
              Transcript_getEnd(t),
              Transcript_getStrand(t),
              Transcript_getBiotype(t),
              Transcript_getAnalysis(t) ? Analysis_getLogicName(Transcript_getAnalysis(t)) : " NO ANALYSIS!");
      if (withSupport) {
        Vector *support = Transcript_getAllSupportingFeatures(t);
        if (support) {
          for (k=0; k<Vector_getNumElement(support); k++) {
            BaseAlignFeature *baf = Vector_getElementAt(support, k);
            fprintf(fp,"   support %s coords: %ld %ld %d evalue %f hCoverage %f\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf), BaseAlignFeature_getpValue(baf), BaseAlignFeature_gethCoverage(baf));
          }
        } else {
          fprintf(fp,"   no transcript support\n");
        }
        Vector *intronSupport = Transcript_getAllIntronSupportingEvidence(t);
        if (intronSupport) {
          for (k=0; k<Vector_getNumElement(intronSupport); k++) {
            IntronSupportingEvidence *ise = Vector_getElementAt(intronSupport, k);
            fprintf(fp,"   intron support %s type %s coords: %ld %ld %d\n",
                    IntronSupportingEvidence_getHitName(ise),
                    Analysis_getLogicName(IntronSupportingEvidence_getAnalysis(ise)),
                    IntronSupportingEvidence_getStart(ise),
                    IntronSupportingEvidence_getEnd(ise),
                    IntronSupportingEvidence_getStrand(ise));
          }
        } else {
          fprintf(fp,"   no intron support\n");
        }
      }

      for (k=0;k<Transcript_getExonCount(t);k++) {
        Exon *e = Transcript_getExonAt(t,k);
        fprintf(fp,"  exon %s (%p) coords: %ld %ld %d\n",Exon_getStableId(e), e, Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
        if (withSupport) {
          Vector *support = Exon_getAllSupportingFeatures(e);
          int m;
          for (m=0; m<Vector_getNumElement(support); m++) {
            BaseAlignFeature *baf = Vector_getElementAt(support, m);
            fprintf(fp,"   support %s coords: %ld %ld %d evalue %f hCoverage %f\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf), BaseAlignFeature_getpValue(baf), BaseAlignFeature_gethCoverage(baf));
          }
        }
      }
      Translation *tln = Transcript_getTranslation(t);
      if (tln) {

        fprintf(fp," translation id: %s %s %d %s %d\n",Translation_getStableId(tln),
                Exon_getStableId(Translation_getStartExon(tln)), Translation_getStart(tln),
                Exon_getStableId(Translation_getEndExon(tln)), Translation_getEnd(tln));
        char *tSeq = Transcript_translate(t);
        fprintf(fp," translation: %s\n",tSeq);
        free(tSeq);
        Vector *tlnAttribs = Translation_getAllAttributes(tln, NULL);
        if (Vector_getNumElement(tlnAttribs)) {
          fprintf(fp, " translation attributes:\n");
          int n;
          for (n=0; n<Vector_getNumElement(tlnAttribs); n++) {
            Attribute *attrib = Vector_getElementAt(tlnAttribs, n);
            fprintf(fp, "  code %s name %s desc %s value %s\n",
                    Attribute_getCode(attrib),
                    Attribute_getName(attrib),
                    Attribute_getDescription(attrib),
                    Attribute_getValue(attrib));
          }
        }
      }
    }
  }
  return failed;
}

// First arg should eventually be RunnableDB type I think - need a method to get the configvar -> setting function hash
void Utilities_parseConfig(RoughTranscripts *rsg, config_setting_t *cfgBlock, char *label, int ignoreThrow) {
  char *DEFAULT_ENTRY_KEY = "DEFAULT";

  if (label == NULL) {
    fprintf(stderr, "Can't parse the config hash for object if we are give no label\n");
    exit(EX_DATAERR);
  }

  StringHash *keyCheckHash = StringHash_new(STRINGHASH_SMALL);

  int count = config_setting_length(cfgBlock);

  int i;
  char ucKey[2048];
  for(i = 0; i < count; ++i) {
    config_setting_t *setting = config_setting_get_elem(cfgBlock, i);

    strcpy(ucKey, config_setting_name(setting));
    StrUtil_strupr(ucKey);

    if (StringHash_contains(keyCheckHash, ucKey)) {
      fprintf(stderr, "You have two entries in your config with the same name (ignoring case) for name %s\n", ucKey);
      exit(EX_DATAERR);
    }
    StringHash_add(keyCheckHash, ucKey, setting);
  }

  config_setting_t *defaultSection = StringHash_getValue(keyCheckHash, DEFAULT_ENTRY_KEY);

  if (defaultSection == NULL) {
    fprintf(stderr, "You must define a %s entry in your config", DEFAULT_ENTRY_KEY);
    exit(EX_DATAERR);
  }

  // the following will fail if there are config variables that
  // do not have a corresponding method here

  count = config_setting_length(defaultSection);

  for(i = 0; i < count; ++i) {
    config_setting_t *setting = config_setting_get_elem(defaultSection, i);

    if (!StringHash_contains(rsg->funcHash, config_setting_name(setting))) {
      fprintf(stderr, "Error: no method defined in Utilities for config variable '%s'\n", config_setting_name(setting));
      exit(EX_DATAERR);
    }

    SetFuncData *setFuncData = StringHash_getValue(rsg->funcHash, config_setting_name(setting));

    ConfigConverter_wrapSetCall(rsg, setFuncData, setting);
  }

  //#########################################################
  // read values of config variables for this logic name into
  // instance variable, set by method
  //#########################################################
  char ucLogic[2048];
  strcpy(ucLogic, label);
  StrUtil_strupr(ucLogic);

  config_setting_t *analSection = StringHash_getValue(keyCheckHash, ucLogic);

  StringHash_free(keyCheckHash, NULL);

  if (analSection != NULL) {
    // entry contains more specific values for the variables
    int count = config_setting_length(analSection);

    int i;
    for(i = 0; i < count; ++i) {
      config_setting_t *setting = config_setting_get_elem(analSection, i);

      if (!StringHash_contains(rsg->funcHash, config_setting_name(setting))) {
        fprintf(stderr, "Error: no method defined in Utilities for config variable '%s'\n", config_setting_name(setting));
        exit(EX_DATAERR);
      }

      SetFuncData *setFuncData = StringHash_getValue(rsg->funcHash, config_setting_name(setting));

      ConfigConverter_wrapSetCall(rsg, setFuncData, setting);
    }

  } else {
    if (ignoreThrow == 1 ) {
      fprintf(stderr,"Warning: Your logic_name %s doesn't appear in your config file hash - using default settings\n",  ucLogic);
    } else {
      fprintf(stderr,"Error: Your logic_name %s doesn't appear in your config file hash - using default settings\n",  ucLogic);
      exit(EX_DATAERR);
    }
  }
}

char *ConfigConverter_typeCodeToString(int code) {
  switch (code) {
    case CONFIG_TYPE_INT:
      return "INT";
      break;

    case CONFIG_TYPE_INT64:
      return "INT64";
      break;

    case CONFIG_TYPE_FLOAT:
      return "FLOAT";
      break;

    case CONFIG_TYPE_STRING:
      return "STRING";
      break;

    case CONFIG_TYPE_BOOL:
      return "BOOL";
      break;

    case CONFIG_TYPE_ARRAY:
      return "ARRAY";
      break;

    case CONFIG_TYPE_LIST:
      return "LIST";
      break;

    case CONFIG_TYPE_GROUP:
      return "GROUP";
      break;

    default:
      fprintf(stderr,"Error: Unknown type in ConfigConverter_typeCodeToString - code = %d\n", code);
      exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapSetCall(RoughTranscripts *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  switch (config_setting_type(setting)) {
    case CONFIG_TYPE_INT:
      if (setFuncData->type != CONFIG_TYPE_INT) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %d), type is %s expected type is %s (config type code)\n",
                config_setting_name(setting), config_setting_get_int(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setIntValue(rsg, config_setting_get_int(setting));
      break;

    case CONFIG_TYPE_INT64:
      if (setFuncData->type != CONFIG_TYPE_INT64) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = "IDFMTSTR"), type is %s expected type is %s (config type code)\n",
                config_setting_name(setting), config_setting_get_int64(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setInt64Value(rsg, config_setting_get_int64(setting));
      break;

    case CONFIG_TYPE_FLOAT:
      if (setFuncData->type != CONFIG_TYPE_FLOAT) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %f), type is %s expected type is %s (config type code)\n",
                config_setting_name(setting), config_setting_get_float(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setFloatValue(rsg, config_setting_get_float(setting));
      break;

    case CONFIG_TYPE_STRING:
      if (setFuncData->type != CONFIG_TYPE_STRING) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %s), type is %s expected type is %s (config type code)\n",
                config_setting_name(setting), config_setting_get_string(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setStringValue(rsg, config_setting_get_string(setting));
      break;

    case CONFIG_TYPE_BOOL:
      if (setFuncData->type != CONFIG_TYPE_BOOL) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %d), type is %s expected type is %s (config type code)\n",
                config_setting_name(setting), config_setting_get_bool(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setBoolValue(rsg, config_setting_get_bool(setting));
      break;

    case CONFIG_TYPE_ARRAY:
      if (setFuncData->type != CONFIG_TYPE_ARRAY) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n",
                config_setting_name(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapArraySetCall(rsg, setFuncData, setting);
      break;

    case CONFIG_TYPE_LIST:
      if (setFuncData->type != CONFIG_TYPE_LIST) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n",
                config_setting_name(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapListSetCall(rsg, setFuncData, setting);
      break;

    case CONFIG_TYPE_GROUP:
      if (setFuncData->type != CONFIG_TYPE_GROUP) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n",
                config_setting_name(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)),
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapGroupSetCall(rsg, setFuncData, setting);
      break;

    default:
      fprintf(stderr,"Error: Unknown type in config type switch for %s, type is %d expected type is %d (config type code)\n",
              config_setting_name(setting), config_setting_type(setting), setFuncData->type);
      exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapGroupSetCall(RoughTranscripts *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  fprintf(stderr,"NIY: Don't know how to handle Group config with name %s\n", config_setting_name(setting));
  exit(EX_DATAERR);
}

void ConfigConverter_wrapListSetCall(RoughTranscripts *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  if (setFuncData->subFunc) {
    setFuncData->subFunc(rsg, setting);
  } else {
    fprintf(stderr,"Don't know how to handle List config with name %s\n", config_setting_name(setting));
    exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapArraySetCall(RoughTranscripts *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  Vector *results = Vector_new();
  if (setFuncData->subType == CONFIG_TYPE_STRING) {
    char *tmp;
    int i;
    for (i=0; i<config_setting_length(setting); i++) {
      const char *val = config_setting_get_string_elem(setting, i);
      Vector_addElement(results, StrUtil_copyString(&tmp, (char*)val, 0));
    }
  } else if (setFuncData->subType == CONFIG_TYPE_FLOAT) {
    int i;
    for (i=0; i<config_setting_length(setting); i++) {
      double val = config_setting_get_float_elem(setting, i);
      Vector_addElement(results, double_new(val));
    }
  } else {
    fprintf(stderr,"Array type other than string not implemented yet\n");
    exit(EX_DATAERR);
  }
  setFuncData->setFunc.setVectorValue(rsg, results);
}


void RoughTranscripts_setDatabaseConfig(RoughTranscripts *rsg, config_setting_t *setting) {
  rsg->databaseConfig = setting;
}

config_setting_t *RoughTranscripts_getDatabaseConfig(RoughTranscripts *rsg) {
  return rsg->databaseConfig;
}

DBAdaptor *BaseGeneBuild_getDbAdaptor(RoughTranscripts *rsg, char *alias, int isNonStandard, int dontUseDnaDb) {
  DBAdaptor *db = NULL;
  int tryToAttachDnaDb = 0;
  int verbosity = RoughTranscripts_getVerbosity(rsg);

  config_setting_t *setting = RoughTranscripts_getDatabaseConfig(rsg);

  config_setting_t *databases = config_setting_get_member(setting, "DATABASES");

  const char *DNA_DBNAME;
  config_setting_lookup_string(setting, "DNA_DBNAME", &DNA_DBNAME);

  if (rsg->adaptorAliasHash == NULL) {
    fprintf(stderr, "Error: adaptorAliasHash is NULL - should have been allocated by now - bye\n");
    exit(EX_SOFTWARE);
  }

  StringHash *hash = rsg->adaptorAliasHash;

  if (!StringHash_contains(hash, alias)) { // if we don't already have an entry for this ...
    config_setting_t *section = config_setting_get_member(databases, alias);
    if (section != NULL) {

      // check if we got all arguments
      const char *user;
      const char *host;
      const char *dbName;
      int   port;
      if (!(config_setting_lookup_string(section, "user", &user) &&
            config_setting_lookup_string(section, "dbname", &dbName) &&
            config_setting_lookup_string(section, "host", &host) &&
            config_setting_lookup_int(section, "port", &port))) {
        fprintf(stderr,"Error: Missing at least one required arg (user, dbname, host or port) in DATABASES config block for %s\n", alias);
        exit(EX_DATAERR);
      }

        const char *pass;
        config_setting_lookup_string(section, "pass", &pass);
        db = DBAdaptor_new((char *)host, (char *)user, (char *)pass, (char *)dbName, port, NULL);
        // it's a core db so try to attach dna_db
        tryToAttachDnaDb = 1;

      // this bit is attaching a dna_db .
      if (DNA_DBNAME != NULL && strcmp(alias, DNA_DBNAME) && tryToAttachDnaDb) {

        if (dontUseDnaDb) {
          if (verbosity > 0) fprintf(stderr,"\nNot attaching a DNA_DB to %s\n.", alias);
          // if two different species are considered, the
          // not_use_dna_database is set to 1 to avoid adding the second species to
          // the first one

        } else {
          // there's a little danger if we have multiple diffeent
          // species in our "Databases.pm" file. We need to avoid that the wrong
          // dna db is attached, ie a mouse core with a human dna db.

          if (verbosity > 0) fprintf(stderr,"\nAttaching DNA_DB %s to %s...\n", DNA_DBNAME, alias);
          if (DNA_DBNAME[0] == '\0') {
            fprintf(stderr, "You're using an empty string as dna_dbname in your Databases config file\n");
            exit(EX_DATAERR);
          }
          DBAdaptor *dnaDb = BaseGeneBuild_getDbAdaptor(rsg, (char *)DNA_DBNAME, 0, 0);

          // try to get default asm+ species name for OTHER db - does not work
          // for comapra database
          MetaContainer *coreMC = DBAdaptor_getMetaContainer(db);

          char *coreDbAsm     = MetaContainer_getDefaultAssembly(coreMC);

          // get the same for dna-db
          MetaContainer *dnaMC = DBAdaptor_getMetaContainer(dnaDb);

          char *dnaDbAsm     = MetaContainer_getDefaultAssembly(dnaMC);

          int dbsAreCompatible = 1;

          if (strcmp(coreDbAsm, dnaDbAsm)) { // assemblies differ
            fprintf(stderr, "You're trying to add a DNA_DB with assembly %s to "
                            "a core/cdna/otherfeatures DB with assembly %s ...\n\t"
                            "that's incompatbile. I will not add dna_database "
                            "%s to core %s\n", dnaDbAsm, coreDbAsm, DBConnection_getDbName(dnaDb->dbc), DBConnection_getDbName(db->dbc));

            dbsAreCompatible = 0;
          }

          free(coreDbAsm);
          free(dnaDbAsm);

          if (dbsAreCompatible) {
            DBAdaptor_setDNADBAdaptor(db, dnaDb);
            fprintf(stderr,"\nAttaching DNA_DB %s to %s\n", DBConnection_getDbName(dnaDb->dbc), DBConnection_getDbName(db->dbc));
          }
        }
      } else {
        if ( !strcmp(alias, DNA_DBNAME)) {
          fprintf(stderr, "\nNot attaching DNA_DB to %s which has DNA_DBNAME...\n", alias);
        } else {
          fprintf(stderr, "You haven't defined a DNA_DBNAME in your Databases config file\n");
        }
      }
    } else {
      fprintf(stderr, "No entry in Databases config file hash for %s\n", alias);
      exit(EX_DATAERR);
    }

    StringHash_add(hash, alias, db);
  } else {
    db = StringHash_getValue(hash, alias);
  }

  return db;
}

void RunnableDB_readDatabaseConfig(RoughTranscripts *rsg, char *configFile) {
  config_t cfg;

  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if (!config_read_file(&cfg, configFile)) {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EX_NOINPUT);
  }

  config_setting_t *cfgBlock = config_lookup(&cfg, "Config");
  if (cfgBlock == NULL) {
    fprintf(stderr,"Missing config block 'Config'\n");
  }

  RoughTranscripts_setDatabaseConfig(rsg, cfgBlock);
}


DAFCluster *DAFCluster_new(DNAAlignFeature *daf) {
  DAFCluster *cluster;
  if ((cluster = (DAFCluster *)calloc(1,sizeof(DAFCluster))) == NULL) {
    fprintf(stderr,"Failed allocating DAFCluster\n");
    exit(EXIT_MEMORY);
  }
  cluster->feature = daf;
  cluster->links = Vector_new();
  cluster->index = -1;
  return cluster;
}

void DAFCluster_setStart(DAFCluster *cluster, int start) {
  DNAAlignFeature_setStart(cluster->feature, start);
}

void DAFCluster_setEnd(DAFCluster *cluster, int end) {
  DNAAlignFeature_setEnd(cluster->feature, end);
}

long DAFCluster_getStart(DAFCluster *cluster) {
  return DNAAlignFeature_getStart(cluster->feature);
}

long DAFCluster_getEnd(DAFCluster *cluster) {
  return DNAAlignFeature_getEnd(cluster->feature);
}

void DAFCluster_incrementHCoverage(DAFCluster *cluster) {
  DNAAlignFeature_sethCoverage(cluster->feature, DNAAlignFeature_gethCoverage(cluster->feature)+1);
}

void DAFCluster_incrementScore(DAFCluster *cluster) {
  DNAAlignFeature_setScore(cluster->feature, DNAAlignFeature_getScore(cluster->feature)+1);
}

void DAFCluster_free(DAFCluster *cluster) {
  Vector_free(cluster->links);
  free(cluster);
}
SetFuncData *SetFuncData_new(void *func, int type) {
  SetFuncData *sfd;

  if ((sfd = (SetFuncData *)calloc(1,sizeof(SetFuncData))) == NULL) {
    fprintf(stderr,"Failed allocating SetFuncData\n");
    exit(EXIT_MEMORY);
  }
  sfd->setFunc.setIntValue = func;
  sfd->type = type;

  return sfd;
}


void SetFuncData_free(SetFuncData *sfd) {
  free(sfd);
}

