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

#ifndef __REFINESOLEXAGENES_H__
#define __REFINESOLEXAGENES_H__

#include <config.h> //Generated by configure

#include "ensc/DBAdaptor.h"
#include "ensc/Vector.h"
#include "ensc/Slice.h"
#include "ensc/StringHash.h"
#include "ensc/SliceAdaptor.h"
#include "ensc/Analysis.h"
#include "ensc/Transcript.h"

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "libconfig.h"

#define RSGEXON_RETAINED 1<<1
#define RSGEXON_EXTRA    1<<2

#define RSGINTRON_NONCANON    1<<3

typedef enum CigarBlockTypeEnum {
  CB_NONE,
  CB_MATCH,
  CB_INTRON,
  CB_DELETION
} CigarBlockType;

typedef struct CigarBlockStruct {
  CigarBlockType type;
  long  start;    // start pos on reference, 1 based I think
  long  end;      // end pos on reference, 1 based I think
} CigarBlock;

typedef struct ExtraExonDataStruct {
  int nCoord;
  long *coords;
  int score;
} ExtraExonData;

typedef struct RefineSolexaGenesStruct {
  char *badModelsType;
  char *bestScoreType;
  char *inputId;
  char *intronDb;
  char *modelDb;
  char *modelLogicName;
  char *outputDb;
  char *otherIsoformsType;
  char *singleExonModelType;
  char *typePrefix;

  int dryRun;
  int max3PrimeExons;
  int max5PrimeExons;
  int max3PrimeLength;
  int max5PrimeLength;
  int maxIntronSize;
  int maxNum;
  int maxRecursions;
  int minIntronSize;
  int minSingleExonLength;
  int otherNum;
  int recursiveLimit;
  int strictInternalSpliceSites;
  int strictInternalEndSpliceSites;
  int trimUtr;
  int verbosity;
  int threads;
  int ucsc_naming;
  int writeIntrons;

  SliceAdaptor *geneSliceAdaptor;
  SliceAdaptor *intronSliceAdaptor;

  Slice *chrSlice;

  StringHash *extraExons;
  StringHash *funcHash;

  char ** extraExonsKeys;

  ExtraExonData *extraExonsValues;

  Vector *intronBamFiles;
  Vector *intronFeatures;
  Vector *logicNames;
  Vector *output;
  Vector *prelimGenes;

  Vector *consLims;
  Vector *nonConsLims;
  double restartNonConsLim;


  int filterOnOverlapThreshold;
  int isOneThreshold;
  double minSingleExonCDSPercLength;
  double rejectIntronCutoff;
  double retainedIntronPenalty;
  double consLim;
  double nonConsLim;

  StringHash *adaptorAliasHash;

  Analysis *analysis;
  Analysis *analysisIntrons;
  Analysis *analysisISE;

  DBAdaptor *db;

  long longestIntronLength;

  config_setting_t *databaseConfig;
} RefineSolexaGenes;

typedef struct ModelClusterStruct {
  long start;
  long end;
  int strand;
  Vector *models; // A vector of Model objects
  Vector *finalModels; // A vector of Gene objects
} ModelCluster;


typedef struct IntronBamConfigStruct {
  int depth;
  int mixedBam;
  Vector *groupNames;
  char *fileName;
} IntronBamConfig;

typedef struct ORFRangeStruct {
  long length;
  long start;
  long end;
} ORFRange;



RefineSolexaGenes *RefineSolexaGenes_new(char *configFile, char *logicName);
DBAdaptor *RefineSolexaGenes_getDbAdaptor(RefineSolexaGenes *rsg, char *alias);
void RefineSolexaGenes_fetchInput(RefineSolexaGenes *rsg);
void  RefineSolexaGenes_run(RefineSolexaGenes *rsg);
void RefineSolexaGenes_refineGenes(RefineSolexaGenes *rsg);
Analysis *RefineSolexaGenes_createAnalysisObject(RefineSolexaGenes *rsg, char *logicName);
Vector *RefineSolexaGenes_reclusterModels(RefineSolexaGenes *rsg, Vector *clusters, Vector **retNewClusters);
ModelCluster *RefineSolexaGenes_recalculateCluster(RefineSolexaGenes *rsg, Vector *genes);
void RefineSolexaGenes_filterModels(RefineSolexaGenes *rsg, Vector *clusters);
Vector *RefineSolexaGenes_makeModels(RefineSolexaGenes *rsg, StringHash *paths, int strand, Vector *exons, Gene *gene, StringHash *intronHash);
Transcript *RefineSolexaGenes_modifyTranscript(RefineSolexaGenes *rsg, Transcript *tran, Vector *exons);
void RefineSolexaGenes_writeOutput(RefineSolexaGenes *rsg);
int RefineSolexaGenes_processTree(RefineSolexaGenes *rsg, StringHash *hashref, char *index, char *soFar, StringHash *paths, int *recLev);
StringHash *RefineSolexaGenes_processPaths(RefineSolexaGenes *rsg, Vector *exons, Vector *exonIntron, StringHash *intronExon, int strict, int *giveUpFlag);
Vector *RefineSolexaGenes_makeModelClusters(RefineSolexaGenes *rsg, Vector *models, int strand);
Vector *RefineSolexaGenes_mergeExons(RefineSolexaGenes *rsg, Gene *gene, int strand);
Exon *RefineSolexaGenes_binSearchForOverlap(RefineSolexaGenes *rsg, Vector *exons, int pos);
void RefineSolexaGenes_bamToIntronFeatures(RefineSolexaGenes *rsg, IntronBamConfig *intronBamConf, htsFile *sam, bam_hdr_t *header, hts_idx_t *idx, int ref, int begRange, int endRange);
int RefineSolexaGenes_getUngappedFeatures(RefineSolexaGenes *rsg, bam_hdr_t *header, bam1_t *b, CigarBlock **ugfs);
void RefineSolexaGenes_dnaToIntronFeatures(RefineSolexaGenes *rsg, long start, long end);
Vector *RefineSolexaGenes_fetchIntronFeatures(RefineSolexaGenes *rsg, long start, long end, long *offsetP);
Exon *RefineSolexaGenes_makeExon(RefineSolexaGenes *rsg, long start, long end, double score, char *diplayId);
int RefineSolexaGenes_pruneUTR(RefineSolexaGenes *rsg, Gene *gene);
void RefineSolexaGenes_setRecursiveLimit(RefineSolexaGenes *rsg, int limit);
int RefineSolexaGenes_getRecursiveLimit(RefineSolexaGenes *rsg);
SliceAdaptor *RefineSolexaGenes_getGeneSliceAdaptor(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setGeneSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa);
SliceAdaptor *RefineSolexaGenes_getIntronSliceAdaptor(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa);
Vector *RefineSolexaGenes_getPrelimGenes(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setPrelimGenes(RefineSolexaGenes *rsg, Vector *genes);
Slice *RefineSolexaGenes_getChrSlice(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setChrSlice(RefineSolexaGenes *rsg, Slice *slice);
Vector *RefineSolexaGenes_getIntronFeatures(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronFeatures(RefineSolexaGenes *rsg, Vector *features);
StringHash *RefineSolexaGenes_getExtraExons(RefineSolexaGenes *rsg);
char **RefineSolexaGenes_getExtraExonsKeys(RefineSolexaGenes *rsg);
ExtraExonData **RefineSolexaGenes_getExtraExonsValues(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setExtraExons(RefineSolexaGenes *rsg, StringHash *extraExons);
void RefineSolexaGenes_setIntronDb(RefineSolexaGenes *rsg, char *intronDb);
char *RefineSolexaGenes_getIntronDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOutputDb(RefineSolexaGenes *rsg, char *outputDb);
char *RefineSolexaGenes_getOutputDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setModelDb(RefineSolexaGenes *rsg, char *modelDb);
char *RefineSolexaGenes_getModelDb(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setLogicNames(RefineSolexaGenes *rsg, Vector *logicNames);
Vector *RefineSolexaGenes_getLogicNames(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setRetainedIntronPenalty(RefineSolexaGenes *rsg, double retainedIntronPenalty);
double RefineSolexaGenes_getRetainedIntronPenalty(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinIntronSize(RefineSolexaGenes *rsg, int minIntronSize);
int RefineSolexaGenes_getMinIntronSize(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxIntronSize(RefineSolexaGenes *rsg, int maxIntronSize);
int RefineSolexaGenes_getMaxIntronSize(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setBestScoreType(RefineSolexaGenes *rsg, char *bestScoreType);
char *RefineSolexaGenes_getBestScoreType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOtherNum(RefineSolexaGenes *rsg, int otherNum);
int RefineSolexaGenes_getOtherNum(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setOtherIsoformsType(RefineSolexaGenes *rsg, char *otherIsoformsType);
char *RefineSolexaGenes_getOtherIsoformsType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setModelLogicName(RefineSolexaGenes *rsg, char *modelLN);
char *RefineSolexaGenes_getModelLogicName(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setBadModelsType(RefineSolexaGenes *rsg, char *badModelsType);
char *RefineSolexaGenes_getBadModelsType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxNum(RefineSolexaGenes *rsg, int maxNum);
int RefineSolexaGenes_getMaxNum(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMaxRecursions(RefineSolexaGenes *rsg, int maxRecursions);
int RefineSolexaGenes_getMaxRecursions(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinSingleExonLength(RefineSolexaGenes *rsg, int minSingleExonLength);
int RefineSolexaGenes_getMinSingleExonLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMinSingleExonCDSPercLength(RefineSolexaGenes *rsg, double minSingleExonCDSPercLength);
double RefineSolexaGenes_getMinSingleExonCDSPercLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setSingleExonModelType(RefineSolexaGenes *rsg, char *singleExonModelType);
char *RefineSolexaGenes_getSingleExonModelType(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setStrictInternalSpliceSites(RefineSolexaGenes *rsg, int strictInternalSpliceSites);
int RefineSolexaGenes_strictInternalSpliceSites(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setStrictInternalEndSpliceSites(RefineSolexaGenes *rsg, int strictInternalEndSpliceSites);
int RefineSolexaGenes_strictInternalEndSpliceSites(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIntronBamFiles(RefineSolexaGenes *rsg, Vector *intronBamFiles);
Vector *RefineSolexaGenes_getIntronBamFiles(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setWriteIntrons(RefineSolexaGenes *rsg, int writeIntrons);
int RefineSolexaGenes_writeIntrons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setTrimUTR(RefineSolexaGenes *rsg, int trimUtr);
int RefineSolexaGenes_trimUTR(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax3PrimeExons(RefineSolexaGenes *rsg, int max3PrimeExons);
int RefineSolexaGenes_getMax3PrimeExons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax3PrimeLength(RefineSolexaGenes *rsg, int max3PrimeLength);
int RefineSolexaGenes_getMax3PrimeLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax5PrimeExons(RefineSolexaGenes *rsg, int max5PrimeExons);
int RefineSolexaGenes_getMax5PrimeExons(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setMax5PrimeLength(RefineSolexaGenes *rsg, int max5PrimeLength);
int RefineSolexaGenes_getMax5PrimeLength(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setFilterOnOverlapThreshold(RefineSolexaGenes *rsg, int filterOnOverlapThreshold);
int RefineSolexaGenes_getFilterOnOverlapThreshold(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setIsOneThreshold(RefineSolexaGenes *rsg, int isOneThreshold);
int RefineSolexaGenes_getIsOneThreshold(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setRejectIntronCutoff(RefineSolexaGenes *rsg, double rejectIntronCutoff);
double RefineSolexaGenes_getRejectIntronCutoff(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setLongestIntronLength(RefineSolexaGenes *rsg, long maxLength);
long RefineSolexaGenes_getLongestIntronLength(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setConsLim(RefineSolexaGenes *rsg, double consLim);
double RefineSolexaGenes_getConsLim(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setNonConsLim(RefineSolexaGenes *rsg, double nonConsLim);
double RefineSolexaGenes_getNonConsLim(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setDatabaseConfig(RefineSolexaGenes *rsg, config_setting_t *setting);
config_setting_t *RefineSolexaGenes_getDatabaseConfig(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setConsLims(RefineSolexaGenes *rsg, Vector *consLims);
Vector *RefineSolexaGenes_getConsLims(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setNonConsLims(RefineSolexaGenes *rsg, Vector *nonConsLims);
Vector *RefineSolexaGenes_getNonConsLims(RefineSolexaGenes *rsg);
void RefineSolexaGenes_setRestartNonConsLim(RefineSolexaGenes *rsg, double restartNonConsLim);
double RefineSolexaGenes_getRestartNonConsLim(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setDryRun(RefineSolexaGenes *rsg, int dryRun);
int RefineSolexaGenes_isDryRun(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setVerbosity(RefineSolexaGenes *rsg, int verbosity);
int RefineSolexaGenes_getVerbosity(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setThreads(RefineSolexaGenes *rsg, int threads);
int RefineSolexaGenes_getThreads(RefineSolexaGenes *rsg);

void RefineSolexaGenes_setUcscNaming(RefineSolexaGenes *rsg, int ucsc_naming);
int RefineSolexaGenes_getUcscNaming(RefineSolexaGenes *rsg);

void RefineSolexaGenes_usage(int exit_code);

void RefineSolexaGenes_dumpConfig(RefineSolexaGenes *rsg);

int SeqFeat_lengthCompFunc(const void *a, const void *b);
void updateDuplicatedGeneBiotype(Gene *gene);
int isGeneDuplicated(RefineSolexaGenes *rsg, Vector *indexes, Vector *genes, Gene *gene);
void RefineSolexaGenes_filterGenes(RefineSolexaGenes *rsg);

// To move
  Transcript *TranslationUtils_addORFToTranscript(ORFRange *orf, Transcript *transcript);
  Transcript *TranslationUtils_computeTranslation(Transcript *tran);
  Gene *TranscriptUtils_convertToGene(Transcript *t, Analysis *analysis, char *biotype);
  Exon *ExonUtils_createExon(long start, long end, int phase, int endPhase, int strand, Analysis *analysis, Vector *supportingFeatures, IDType dbId, Slice *slice, char *stableId, int version);
  Exon *ExonUtils_cloneExon(Exon *exon);
  BaseAlignFeature *EvidenceUtils_cloneEvidence(BaseAlignFeature *feature);
  Slice *RefineSolexaGenes_fetchSequence(RefineSolexaGenes *rsg, char *name, DBAdaptor *db, Vector *repeatMaskTypes, int softMask);
  void RefineSolexaGenes_setInputId(RefineSolexaGenes *rsg, char *inputId);
  char *RefineSolexaGenes_getInputId(RefineSolexaGenes *rsg);
  void RefineSolexaGenes_setTypePrefix(RefineSolexaGenes *rsg, char *typePrefix);
  char *RefineSolexaGenes_getTypePrefix(RefineSolexaGenes *rsg);
  Vector *RefineSolexaGenes_getOutput(RefineSolexaGenes *rsg);
  void RefineSolexaGenes_addToOutput(RefineSolexaGenes *rsg, Gene *gene);
  void RefineSolexaGenes_setDb(RefineSolexaGenes *rsg, DBAdaptor *db);
  DBAdaptor *RefineSolexaGenes_getDb(RefineSolexaGenes *rsg);
  Vector *TranslationUtils_generateORFRanges(Transcript *transcript, int requireMet, int minLength, int allowReverse);

  void RefineSolexaGenes_setAnalysis(RefineSolexaGenes *rsg, Analysis *analysis);
  Analysis *RefineSolexaGenes_getAnalysis(RefineSolexaGenes *rsg);

  void RefineSolexaGenes_setIntronsAnalysis(RefineSolexaGenes *rsg, char *intronsLn);
  Analysis *RefineSolexaGenes_getIntronsAnalysis(RefineSolexaGenes *rsg);

  void RefineSolexaGenes_setISEAnalysis(RefineSolexaGenes *rsg, char *iseLn);
  Analysis *RefineSolexaGenes_getISEAnalysis(RefineSolexaGenes *rsg);


  typedef (*IntSetFunc)(RefineSolexaGenes *rsg, int val);
  typedef (*Int64SetFunc)(RefineSolexaGenes *rsg, long val);
  typedef (*FloatSetFunc)(RefineSolexaGenes *rsg, double val);
  typedef (*StringSetFunc)(RefineSolexaGenes *rsg, const char *val);
  typedef (*BoolSetFunc)(RefineSolexaGenes *rsg, int val);
  typedef (*VectorSetFunc)(RefineSolexaGenes *rsg, Vector *val);
  typedef union {
    IntSetFunc setIntValue;
    Int64SetFunc setInt64Value;
    FloatSetFunc setFloatValue;
    StringSetFunc setStringValue;
    BoolSetFunc setBoolValue;
    VectorSetFunc setVectorValue;
  } SetFunc;

  typedef (*SetSubFunc)(RefineSolexaGenes *rsg, config_setting_t *subSetting);
  typedef struct SetFuncDataStruct {
    int type;
    int subType;
    SetSubFunc subFunc;
    SetFunc setFunc;
  } SetFuncData;

  void RefineSolexaGenes_parseIntronBamFilesConfig(RefineSolexaGenes *rsg, config_setting_t *setting);
  void Utilities_parseConfig(RefineSolexaGenes *rsg, config_setting_t *cfgBlock, char *label, int ignoreThrow);
  char *ConfigConverter_typeCodeToString(int code);
  void ConfigConverter_wrapSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting);
  void ConfigConverter_wrapGroupSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting);
  void ConfigConverter_wrapArraySetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting);
  void ConfigConverter_wrapListSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting);
  void RunnableDB_readAndCheckConfig(RefineSolexaGenes *rsg, char *configFile, char *blockName, char *logicName);

  DBAdaptor *BaseGeneBuild_getDbAdaptor(RefineSolexaGenes *rsg, char *alias, int isNonStandard, int dontUseDnaDb);
  void RunnableDB_readDatabaseConfig(RefineSolexaGenes *rsg, char *configFile);
#endif