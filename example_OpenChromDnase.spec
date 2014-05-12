#
# Example specification file for TACO.
#
# Comprehensive prediction of transcription factor dimers in cell-type-specific open chromatin regions.

<Genome>
  FastaFile = hg19/*.fa
  MaskedRegions = coding_hg19.bed
</Genome>

#
# Open chromatin datasets, e.g. DNase-seq peaks
#
# One replicate per line, two fields are required: dataset name and BED filename.
# For each dataset, the union of all corresponding replicates will be taken.

<StronglySpecificDatasets>
  DatasetList = wgEncodeOpenChromDnase_hg19.list

  # Dataset normalization (each replicate separately):
  # make all hypersensitive regions the same size
  RegionSize = 300
  # exclude a hypersensitive region if most of the underlying genomic sequence is masked
  RegionMasking = Majority
  # consider not more than the given number of regions with top signalValue
  RegionCount = 50000
</StronglySpecificDatasets>

#
# Motif database, e.g. TRANSFAC, JASPAR, SwissRegulon -- preferably use only one of these
#

# TRANSFAC
<Motifs>
  Database = TRANSFAC/matrix.dat
  DatabaseSubset = TRANSFAC.vertebrata
  Sensitivity = 0.8
</Motifs>

# JASPAR
#<Motifs>
#  Database = JASPAR/jaspar_CORE/non_redundant/by_tax_group/vertebrates/matrix_only/matrix_only.txt
#  Sensitivity = 0.9
#</Motifs>

# SwissRegulon
#<Motifs>
#  Database = SwissRegulon/weight_matrices
#  Sensitivity = 0.95
#</Motifs>

#
# Various options, do not forget to adjust NumberOfThreads
#

<Options>
  NumberOfThreads = 16

  MinMotifInformationContribution = 6.0
  MaxOverlappingInformationContent = 2.0
  MaxMotifSpacing = 50
  ConsiderOrientationsSeparately = True
  ConsiderMostSignificantComplexOnly = False

  TargetInstancesThreshold = 100
  FoldChangeThreshold = 1.0
  PValueThreshold = 0.05

  DimerMotifFlanks = 5
  ClusteringAcrossDatasets = True
  ClusteringDistanceConstant = 0.0
  ClusteringDistanceMultiplier = 0.15
  ClusteringOverlapThreshold = 0.2

  OutputDetailedStats = All
  OutputDimerMotifs = All
  OutputGenomicLocations = All
  GenomicLocationsMaxSpacingDeviation = 0
  OutputPValueDistribution = True
</Options>
