#
# Example specification file for TACO.
#
# Comprehensive prediction of transcription factor dimers in K562 ChIP-seq peaks.

<Genome>
  FastaFile = hg19/*.fa
  MaskedRegions = coding_hg19.bed
</Genome>

#
# ChIP-seq datasets
#

<WeaklySpecificDatasets>
  DatasetList = K562_ChIP-seq_hg19.list
</WeaklySpecificDatasets>

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

# Additional motifs, e.g. discovered by MEME (in PSPM format)
<Motifs>
  Motif = Broad/CTCF:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeBroadHistoneK562CtcfStdAlnRep0/motif1
  Motif = Broad/CTCF:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeBroadHistoneK562CtcfStdAlnRep0/motif2
  Motif = Broad/CTCF:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeBroadHistoneK562CtcfStdAlnRep0/motif3
  Motif = Broad/CTCF:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeBroadHistoneK562CtcfStdAlnRep0/motif4
  Motif = Broad/CTCF:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeBroadHistoneK562CtcfStdAlnRep0/motif5

  Motif = HudsonAlpha/ATF3:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif1
  Motif = HudsonAlpha/ATF3:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif2
  Motif = HudsonAlpha/ATF3:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif3
  Motif = HudsonAlpha/ATF3:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif4
  Motif = HudsonAlpha/ATF3:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif5

  Motif = HudsonAlpha/BCLAF1_(SC-101388):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Bclaf101388Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/BCLAF1_(SC-101388):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Bclaf101388Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/BCLAF1_(SC-101388):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Bclaf101388Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/BCLAF1_(SC-101388):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Bclaf101388Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/BCLAF1_(SC-101388):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Bclaf101388Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/CTCF_(SC-5916):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif1
  Motif = HudsonAlpha/CTCF_(SC-5916):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif2
  Motif = HudsonAlpha/CTCF_(SC-5916):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif3
  Motif = HudsonAlpha/CTCF_(SC-5916):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif4
  Motif = HudsonAlpha/CTCF_(SC-5916):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif5

  Motif = HudsonAlpha/CTCFL_(SC-98982):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ctcflsc98982V0416101AlnRep0/motif1
  Motif = HudsonAlpha/CTCFL_(SC-98982):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ctcflsc98982V0416101AlnRep0/motif2
  Motif = HudsonAlpha/CTCFL_(SC-98982):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ctcflsc98982V0416101AlnRep0/motif3
  Motif = HudsonAlpha/CTCFL_(SC-98982):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ctcflsc98982V0416101AlnRep0/motif4
  Motif = HudsonAlpha/CTCFL_(SC-98982):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ctcflsc98982V0416101AlnRep0/motif5

  Motif = HudsonAlpha/E2F6:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562E2f6sc22823V0416102AlnRep0/motif1
  Motif = HudsonAlpha/E2F6:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562E2f6sc22823V0416102AlnRep0/motif2
  Motif = HudsonAlpha/E2F6:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562E2f6sc22823V0416102AlnRep0/motif3
  Motif = HudsonAlpha/E2F6:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562E2f6sc22823V0416102AlnRep0/motif4
  Motif = HudsonAlpha/E2F6:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562E2f6sc22823V0416102AlnRep0/motif5

  Motif = HudsonAlpha/Egr-1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Egr1V0416101AlnRep0/motif1
  Motif = HudsonAlpha/Egr-1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Egr1V0416101AlnRep0/motif2
  Motif = HudsonAlpha/Egr-1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Egr1V0416101AlnRep0/motif3
  Motif = HudsonAlpha/Egr-1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Egr1V0416101AlnRep0/motif4
  Motif = HudsonAlpha/Egr-1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Egr1V0416101AlnRep0/motif5

  Motif = HudsonAlpha/ELF1_(SC-631):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Elf1sc631V0416102AlnRep0/motif1
  Motif = HudsonAlpha/ELF1_(SC-631):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Elf1sc631V0416102AlnRep0/motif2
  Motif = HudsonAlpha/ELF1_(SC-631):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Elf1sc631V0416102AlnRep0/motif3
  Motif = HudsonAlpha/ELF1_(SC-631):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Elf1sc631V0416102AlnRep0/motif4
  Motif = HudsonAlpha/ELF1_(SC-631):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Elf1sc631V0416102AlnRep0/motif5

  Motif = HudsonAlpha/ETS1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ets1V0416101AlnRep0/motif1
  Motif = HudsonAlpha/ETS1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ets1V0416101AlnRep0/motif2
  Motif = HudsonAlpha/ETS1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ets1V0416101AlnRep0/motif3
  Motif = HudsonAlpha/ETS1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ets1V0416101AlnRep0/motif4
  Motif = HudsonAlpha/ETS1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Ets1V0416101AlnRep0/motif5

  Motif = HudsonAlpha/FOSL1_(SC-183):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Fosl1sc183V0416101AlnRep0/motif1
  Motif = HudsonAlpha/FOSL1_(SC-183):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Fosl1sc183V0416101AlnRep0/motif2
  Motif = HudsonAlpha/FOSL1_(SC-183):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Fosl1sc183V0416101AlnRep0/motif3
  Motif = HudsonAlpha/FOSL1_(SC-183):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Fosl1sc183V0416101AlnRep0/motif4
  Motif = HudsonAlpha/FOSL1_(SC-183):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Fosl1sc183V0416101AlnRep0/motif5

  Motif = HudsonAlpha/GABP:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562GabpV0416101AlnRep0/motif1
  Motif = HudsonAlpha/GABP:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562GabpV0416101AlnRep0/motif2
  Motif = HudsonAlpha/GABP:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562GabpV0416101AlnRep0/motif3
  Motif = HudsonAlpha/GABP:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562GabpV0416101AlnRep0/motif4
  Motif = HudsonAlpha/GABP:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562GabpV0416101AlnRep0/motif5

  Motif = HudsonAlpha/GATA2_(SC-267):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Gata2sc267Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/GATA2_(SC-267):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Gata2sc267Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/GATA2_(SC-267):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Gata2sc267Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/GATA2_(SC-267):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Gata2sc267Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/GATA2_(SC-267):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Gata2sc267Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/HDAC2_(SC-6296):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Hdac2sc6296V0416102AlnRep0/motif1
  Motif = HudsonAlpha/HDAC2_(SC-6296):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Hdac2sc6296V0416102AlnRep0/motif2
  Motif = HudsonAlpha/HDAC2_(SC-6296):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Hdac2sc6296V0416102AlnRep0/motif3
  Motif = HudsonAlpha/HDAC2_(SC-6296):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Hdac2sc6296V0416102AlnRep0/motif4
  Motif = HudsonAlpha/HDAC2_(SC-6296):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Hdac2sc6296V0416102AlnRep0/motif5

  Motif = HudsonAlpha/Max:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562MaxV0416102AlnRep0/motif1
  Motif = HudsonAlpha/Max:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562MaxV0416102AlnRep0/motif2
  Motif = HudsonAlpha/Max:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562MaxV0416102AlnRep0/motif3
  Motif = HudsonAlpha/Max:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562MaxV0416102AlnRep0/motif4
  Motif = HudsonAlpha/Max:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562MaxV0416102AlnRep0/motif5

  Motif = HudsonAlpha/MEF2A:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsGm12878Mef2aPcr1xAlnRep0/motif1
  Motif = HudsonAlpha/MEF2A:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsGm12878Mef2aPcr1xAlnRep0/motif2
  Motif = HudsonAlpha/MEF2A:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsGm12878Mef2aPcr1xAlnRep0/motif3
  Motif = HudsonAlpha/MEF2A:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsGm12878Mef2aPcr1xAlnRep0/motif4
  Motif = HudsonAlpha/MEF2A:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsGm12878Mef2aPcr1xAlnRep0/motif5

  Motif = HudsonAlpha/NRSF:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562NrsfV0416102AlnRep0/motif1
  Motif = HudsonAlpha/NRSF:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562NrsfV0416102AlnRep0/motif2
  Motif = HudsonAlpha/NRSF:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562NrsfV0416102AlnRep0/motif3
  Motif = HudsonAlpha/NRSF:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562NrsfV0416102AlnRep0/motif4
  Motif = HudsonAlpha/NRSF:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562NrsfV0416102AlnRep0/motif5

  Motif = HudsonAlpha/PU.1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/PU.1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/PU.1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/PU.1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/PU.1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/Rad21:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Rad21V0416102AlnRep0/motif1
  Motif = HudsonAlpha/Rad21:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Rad21V0416102AlnRep0/motif2
  Motif = HudsonAlpha/Rad21:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Rad21V0416102AlnRep0/motif3
  Motif = HudsonAlpha/Rad21:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Rad21V0416102AlnRep0/motif4
  Motif = HudsonAlpha/Rad21:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Rad21V0416102AlnRep0/motif5

  Motif = HudsonAlpha/Sin3Ak-20:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sin3ak20V0416101AlnRep0/motif1
  Motif = HudsonAlpha/Sin3Ak-20:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sin3ak20V0416101AlnRep0/motif2
  Motif = HudsonAlpha/Sin3Ak-20:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sin3ak20V0416101AlnRep0/motif3
  Motif = HudsonAlpha/Sin3Ak-20:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sin3ak20V0416101AlnRep0/motif4
  Motif = HudsonAlpha/Sin3Ak-20:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sin3ak20V0416101AlnRep0/motif5

  Motif = HudsonAlpha/SIX5:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Six5Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/SIX5:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Six5Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/SIX5:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Six5Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/SIX5:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Six5Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/SIX5:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Six5Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/SP1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp1Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/SP1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp1Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/SP1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp1Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/SP1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp1Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/SP1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp1Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/SP2_(SC-643):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp2sc643V0416102AlnRep0/motif1
  Motif = HudsonAlpha/SP2_(SC-643):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp2sc643V0416102AlnRep0/motif2
  Motif = HudsonAlpha/SP2_(SC-643):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp2sc643V0416102AlnRep0/motif3
  Motif = HudsonAlpha/SP2_(SC-643):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp2sc643V0416102AlnRep0/motif4
  Motif = HudsonAlpha/SP2_(SC-643):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Sp2sc643V0416102AlnRep0/motif5

  Motif = HudsonAlpha/SRF:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562SrfV0416101AlnRep0/motif1
  Motif = HudsonAlpha/SRF:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562SrfV0416101AlnRep0/motif2
  Motif = HudsonAlpha/SRF:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562SrfV0416101AlnRep0/motif3
  Motif = HudsonAlpha/SRF:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562SrfV0416101AlnRep0/motif4
  Motif = HudsonAlpha/SRF:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562SrfV0416101AlnRep0/motif5

  Motif = HudsonAlpha/TAF1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf1V0416101AlnRep0/motif1
  Motif = HudsonAlpha/TAF1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf1V0416101AlnRep0/motif2
  Motif = HudsonAlpha/TAF1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf1V0416101AlnRep0/motif3
  Motif = HudsonAlpha/TAF1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf1V0416101AlnRep0/motif4
  Motif = HudsonAlpha/TAF1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf1V0416101AlnRep0/motif5

  Motif = HudsonAlpha/TAF7_(SC-101167):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf7sc101167V0416101AlnRep0/motif1
  Motif = HudsonAlpha/TAF7_(SC-101167):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf7sc101167V0416101AlnRep0/motif2
  Motif = HudsonAlpha/TAF7_(SC-101167):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf7sc101167V0416101AlnRep0/motif3
  Motif = HudsonAlpha/TAF7_(SC-101167):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf7sc101167V0416101AlnRep0/motif4
  Motif = HudsonAlpha/TAF7_(SC-101167):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Taf7sc101167V0416101AlnRep0/motif5

  Motif = HudsonAlpha/THAP1_(SC-98174):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Thap1sc98174V0416101AlnRep0/motif1
  Motif = HudsonAlpha/THAP1_(SC-98174):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Thap1sc98174V0416101AlnRep0/motif2
  Motif = HudsonAlpha/THAP1_(SC-98174):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Thap1sc98174V0416101AlnRep0/motif3
  Motif = HudsonAlpha/THAP1_(SC-98174):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Thap1sc98174V0416101AlnRep0/motif4
  Motif = HudsonAlpha/THAP1_(SC-98174):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Thap1sc98174V0416101AlnRep0/motif5

  Motif = HudsonAlpha/USF-1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Usf1V0416101AlnRep0/motif1
  Motif = HudsonAlpha/USF-1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Usf1V0416101AlnRep0/motif2
  Motif = HudsonAlpha/USF-1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Usf1V0416101AlnRep0/motif3
  Motif = HudsonAlpha/USF-1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Usf1V0416101AlnRep0/motif4
  Motif = HudsonAlpha/USF-1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Usf1V0416101AlnRep0/motif5

  Motif = HudsonAlpha/YY1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416102AlnRep0/motif1
  Motif = HudsonAlpha/YY1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416102AlnRep0/motif2
  Motif = HudsonAlpha/YY1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416102AlnRep0/motif3
  Motif = HudsonAlpha/YY1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416102AlnRep0/motif4
  Motif = HudsonAlpha/YY1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416102AlnRep0/motif5

  Motif = HudsonAlpha/YY1_(SC-281):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416101AlnRep0/motif1
  Motif = HudsonAlpha/YY1_(SC-281):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416101AlnRep0/motif2
  Motif = HudsonAlpha/YY1_(SC-281):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416101AlnRep0/motif3
  Motif = HudsonAlpha/YY1_(SC-281):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416101AlnRep0/motif4
  Motif = HudsonAlpha/YY1_(SC-281):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Yy1V0416101AlnRep0/motif5

  Motif = HudsonAlpha/ZBTB33:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb33Pcr1xAlnRep0/motif1
  Motif = HudsonAlpha/ZBTB33:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb33Pcr1xAlnRep0/motif2
  Motif = HudsonAlpha/ZBTB33:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb33Pcr1xAlnRep0/motif3
  Motif = HudsonAlpha/ZBTB33:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb33Pcr1xAlnRep0/motif4
  Motif = HudsonAlpha/ZBTB33:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb33Pcr1xAlnRep0/motif5

  Motif = HudsonAlpha/ZBTB7A_(SC-34508):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb7asc34508V0416101AlnRep0/motif1
  Motif = HudsonAlpha/ZBTB7A_(SC-34508):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb7asc34508V0416101AlnRep0/motif2
  Motif = HudsonAlpha/ZBTB7A_(SC-34508):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb7asc34508V0416101AlnRep0/motif3
  Motif = HudsonAlpha/ZBTB7A_(SC-34508):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb7asc34508V0416101AlnRep0/motif4
  Motif = HudsonAlpha/ZBTB7A_(SC-34508):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsK562Zbtb7asc34508V0416101AlnRep0/motif5

  Motif = Harvard/ATF3:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif1
  Motif = Harvard/ATF3:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif2
  Motif = Harvard/ATF3:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif3
  Motif = Harvard/ATF3:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif4
  Motif = Harvard/ATF3:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Atf3StdAlnRep0/motif5

  Motif = Harvard/BDP1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Bdp1StdAlnRep0/motif1
  Motif = Harvard/BDP1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Bdp1StdAlnRep0/motif2
  Motif = Harvard/BDP1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Bdp1StdAlnRep0/motif3
  Motif = Harvard/BDP1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Bdp1StdAlnRep0/motif4
  Motif = Harvard/BDP1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Bdp1StdAlnRep0/motif5

  Motif = Stanford/BHLHE40_(NB100-1800):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsHepg2Bhlhe40V0416101AlnRep0/motif1
  Motif = Stanford/BHLHE40_(NB100-1800):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsHepg2Bhlhe40V0416101AlnRep0/motif2
  Motif = Stanford/BHLHE40_(NB100-1800):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsHepg2Bhlhe40V0416101AlnRep0/motif3
  Motif = Stanford/BHLHE40_(NB100-1800):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsHepg2Bhlhe40V0416101AlnRep0/motif4
  Motif = Stanford/BHLHE40_(NB100-1800):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeHaibTfbsHepg2Bhlhe40V0416101AlnRep0/motif5

  Motif = Harvard/BRF1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf1StdAlnRep0/motif1
  Motif = Harvard/BRF1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf1StdAlnRep0/motif2
  Motif = Harvard/BRF1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf1StdAlnRep0/motif3
  Motif = Harvard/BRF1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf1StdAlnRep0/motif4
  Motif = Harvard/BRF1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf1StdAlnRep0/motif5

  Motif = Harvard/BRF2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf2StdAlnRep0/motif1
  Motif = Harvard/BRF2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf2StdAlnRep0/motif2
  Motif = Harvard/BRF2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf2StdAlnRep0/motif3
  Motif = Harvard/BRF2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf2StdAlnRep0/motif4
  Motif = Harvard/BRF2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brf2StdAlnRep0/motif5

  Motif = Stanford/Brg1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brg1IggmusAlnRep0/motif1
  Motif = Stanford/Brg1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brg1IggmusAlnRep0/motif2
  Motif = Stanford/Brg1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brg1IggmusAlnRep0/motif3
  Motif = Stanford/Brg1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brg1IggmusAlnRep0/motif4
  Motif = Stanford/Brg1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Brg1IggmusAlnRep0/motif5

  Motif = Harvard/CCNT2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ccnt2StdAlnRep0/motif1
  Motif = Harvard/CCNT2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ccnt2StdAlnRep0/motif2
  Motif = Harvard/CCNT2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ccnt2StdAlnRep0/motif3
  Motif = Harvard/CCNT2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ccnt2StdAlnRep0/motif4
  Motif = Harvard/CCNT2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ccnt2StdAlnRep0/motif5

  Motif = Stanford/CEBPB:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CebpbIggrabAlnRep0/motif1
  Motif = Stanford/CEBPB:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CebpbIggrabAlnRep0/motif2
  Motif = Stanford/CEBPB:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CebpbIggrabAlnRep0/motif3
  Motif = Stanford/CEBPB:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CebpbIggrabAlnRep0/motif4
  Motif = Stanford/CEBPB:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CebpbIggrabAlnRep0/motif5

  Motif = Yale/c-Fos:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CfosStdAlnRep0/motif1
  Motif = Yale/c-Fos:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CfosStdAlnRep0/motif2
  Motif = Yale/c-Fos:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CfosStdAlnRep0/motif3
  Motif = Yale/c-Fos:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CfosStdAlnRep0/motif4
  Motif = Yale/c-Fos:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CfosStdAlnRep0/motif5

  Motif = Stanford/CHD2_(AB68301):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Chd2ab68301IggrabAlnRep0/motif1
  Motif = Stanford/CHD2_(AB68301):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Chd2ab68301IggrabAlnRep0/motif2
  Motif = Stanford/CHD2_(AB68301):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Chd2ab68301IggrabAlnRep0/motif3
  Motif = Stanford/CHD2_(AB68301):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Chd2ab68301IggrabAlnRep0/motif4
  Motif = Stanford/CHD2_(AB68301):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Chd2ab68301IggrabAlnRep0/motif5

  Motif = Yale/c-Jun/IFNa6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfna6hStdAlnRep0/motif1
  Motif = Yale/c-Jun/IFNa6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfna6hStdAlnRep0/motif2
  Motif = Yale/c-Jun/IFNa6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfna6hStdAlnRep0/motif3
  Motif = Yale/c-Jun/IFNa6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfna6hStdAlnRep0/motif4
  Motif = Yale/c-Jun/IFNa6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfna6hStdAlnRep0/motif5

  Motif = Yale/c-Jun/IFNg30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng30StdAlnRep0/motif1
  Motif = Yale/c-Jun/IFNg30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng30StdAlnRep0/motif2
  Motif = Yale/c-Jun/IFNg30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng30StdAlnRep0/motif3
  Motif = Yale/c-Jun/IFNg30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng30StdAlnRep0/motif4
  Motif = Yale/c-Jun/IFNg30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng30StdAlnRep0/motif5

  Motif = Yale/c-Jun/IFNg6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng6hStdAlnRep0/motif1
  Motif = Yale/c-Jun/IFNg6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng6hStdAlnRep0/motif2
  Motif = Yale/c-Jun/IFNg6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng6hStdAlnRep0/motif3
  Motif = Yale/c-Jun/IFNg6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng6hStdAlnRep0/motif4
  Motif = Yale/c-Jun/IFNg6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunIfng6hStdAlnRep0/motif5

  Motif = Yale/c-Jun:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunStdAlnRep0/motif1
  Motif = Yale/c-Jun:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunStdAlnRep0/motif2
  Motif = Yale/c-Jun:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunStdAlnRep0/motif3
  Motif = Yale/c-Jun:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunStdAlnRep0/motif4
  Motif = Yale/c-Jun:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CjunStdAlnRep0/motif5

  Motif = Yale/c-Myc/IFNa30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna30StdAlnRep0/motif1
  Motif = Yale/c-Myc/IFNa30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna30StdAlnRep0/motif2
  Motif = Yale/c-Myc/IFNa30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna30StdAlnRep0/motif3
  Motif = Yale/c-Myc/IFNa30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna30StdAlnRep0/motif4
  Motif = Yale/c-Myc/IFNa30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna30StdAlnRep0/motif5

  Motif = Yale/c-Myc/IFNa6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna6hStdAlnRep0/motif1
  Motif = Yale/c-Myc/IFNa6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna6hStdAlnRep0/motif2
  Motif = Yale/c-Myc/IFNa6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna6hStdAlnRep0/motif3
  Motif = Yale/c-Myc/IFNa6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna6hStdAlnRep0/motif4
  Motif = Yale/c-Myc/IFNa6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfna6hStdAlnRep0/motif5

  Motif = Stanford/c-Myc/IFNg30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng30StdAlnRep0/motif1
  Motif = Stanford/c-Myc/IFNg30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng30StdAlnRep0/motif2
  Motif = Stanford/c-Myc/IFNg30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng30StdAlnRep0/motif3
  Motif = Stanford/c-Myc/IFNg30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng30StdAlnRep0/motif4
  Motif = Stanford/c-Myc/IFNg30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng30StdAlnRep0/motif5

  Motif = Yale/c-Myc/IFNg6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng6hStdAlnRep0/motif1
  Motif = Yale/c-Myc/IFNg6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng6hStdAlnRep0/motif2
  Motif = Yale/c-Myc/IFNg6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng6hStdAlnRep0/motif3
  Motif = Yale/c-Myc/IFNg6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng6hStdAlnRep0/motif4
  Motif = Yale/c-Myc/IFNg6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycIfng6hStdAlnRep0/motif5

  Motif = Yale/c-Myc:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycStdAlnRep0/motif1
  Motif = Yale/c-Myc:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycStdAlnRep0/motif2
  Motif = Yale/c-Myc:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycStdAlnRep0/motif3
  Motif = Yale/c-Myc:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycStdAlnRep0/motif4
  Motif = Yale/c-Myc:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562CmycStdAlnRep0/motif5

  Motif = USC/E2F4:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f4UcdAlnRep0/motif1
  Motif = USC/E2F4:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f4UcdAlnRep0/motif2
  Motif = USC/E2F4:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f4UcdAlnRep0/motif3
  Motif = USC/E2F4:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f4UcdAlnRep0/motif4
  Motif = USC/E2F4:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f4UcdAlnRep0/motif5

  Motif = USC/E2F6:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f6UcdAlnRep0/motif1
  Motif = USC/E2F6:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f6UcdAlnRep0/motif2
  Motif = USC/E2F6:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f6UcdAlnRep0/motif3
  Motif = USC/E2F6:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f6UcdAlnRep0/motif4
  Motif = USC/E2F6:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562E2f6UcdAlnRep0/motif5

  Motif = USC/GATA-1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata1UcdAlnRep0/motif1
  Motif = USC/GATA-1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata1UcdAlnRep0/motif2
  Motif = USC/GATA-1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata1UcdAlnRep0/motif3
  Motif = USC/GATA-1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata1UcdAlnRep0/motif4
  Motif = USC/GATA-1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata1UcdAlnRep0/motif5

  Motif = USC/GATA-2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata2UcdAlnRep0/motif1
  Motif = USC/GATA-2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata2UcdAlnRep0/motif2
  Motif = USC/GATA-2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata2UcdAlnRep0/motif3
  Motif = USC/GATA-2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata2UcdAlnRep0/motif4
  Motif = USC/GATA-2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gata2UcdAlnRep0/motif5

  Motif = Harvard/GTF2B:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2bStdAlnRep0/motif1
  Motif = Harvard/GTF2B:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2bStdAlnRep0/motif2
  Motif = Harvard/GTF2B:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2bStdAlnRep0/motif3
  Motif = Harvard/GTF2B:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2bStdAlnRep0/motif4
  Motif = Harvard/GTF2B:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2bStdAlnRep0/motif5

  Motif = Stanford/GTF2F1_(AB28179):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2f1ab28179IggrabAlnRep0/motif1
  Motif = Stanford/GTF2F1_(AB28179):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2f1ab28179IggrabAlnRep0/motif2
  Motif = Stanford/GTF2F1_(AB28179):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2f1ab28179IggrabAlnRep0/motif3
  Motif = Stanford/GTF2F1_(AB28179):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2f1ab28179IggrabAlnRep0/motif4
  Motif = Stanford/GTF2F1_(AB28179):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Gtf2f1ab28179IggrabAlnRep0/motif5

  Motif = Harvard/HMGN3:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Hmgn3StdAlnRep0/motif1
  Motif = Harvard/HMGN3:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Hmgn3StdAlnRep0/motif2
  Motif = Harvard/HMGN3:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Hmgn3StdAlnRep0/motif3
  Motif = Harvard/HMGN3:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Hmgn3StdAlnRep0/motif4
  Motif = Harvard/HMGN3:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Hmgn3StdAlnRep0/motif5

  Motif = Stanford/Ini1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ini1IggmusAlnRep0/motif1
  Motif = Stanford/Ini1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ini1IggmusAlnRep0/motif2
  Motif = Stanford/Ini1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ini1IggmusAlnRep0/motif3
  Motif = Stanford/Ini1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ini1IggmusAlnRep0/motif4
  Motif = Stanford/Ini1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Ini1IggmusAlnRep0/motif5

  Motif = Stanford/IRF1/IFNa30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifna30StdAlnRep0/motif1
  Motif = Stanford/IRF1/IFNa30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifna30StdAlnRep0/motif2
  Motif = Stanford/IRF1/IFNa30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifna30StdAlnRep0/motif3
  Motif = Stanford/IRF1/IFNa30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifna30StdAlnRep0/motif4
  Motif = Stanford/IRF1/IFNa30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifna30StdAlnRep0/motif5

  Motif = Stanford/IRF1/IFNg6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifng6hStdAlnRep0/motif1
  Motif = Stanford/IRF1/IFNg6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifng6hStdAlnRep0/motif2
  Motif = Stanford/IRF1/IFNg6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifng6hStdAlnRep0/motif3
  Motif = Stanford/IRF1/IFNg6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifng6hStdAlnRep0/motif4
  Motif = Stanford/IRF1/IFNg6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Irf1Ifng6hStdAlnRep0/motif5

  Motif = Stanford/JunD:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562JundStdAlnRep0/motif1
  Motif = Stanford/JunD:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562JundStdAlnRep0/motif2
  Motif = Stanford/JunD:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562JundStdAlnRep0/motif3
  Motif = Stanford/JunD:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562JundStdAlnRep0/motif4
  Motif = Stanford/JunD:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562JundStdAlnRep0/motif5

  Motif = Stanford/MafK_(ab50322):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mafkab50322IggrabAlnRep0/motif1
  Motif = Stanford/MafK_(ab50322):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mafkab50322IggrabAlnRep0/motif2
  Motif = Stanford/MafK_(ab50322):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mafkab50322IggrabAlnRep0/motif3
  Motif = Stanford/MafK_(ab50322):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mafkab50322IggrabAlnRep0/motif4
  Motif = Stanford/MafK_(ab50322):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mafkab50322IggrabAlnRep0/motif5

  Motif = Yale/Max:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562MaxStdAlnRep0/motif1
  Motif = Yale/Max:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562MaxStdAlnRep0/motif2
  Motif = Yale/Max:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562MaxStdAlnRep0/motif3
  Motif = Yale/Max:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562MaxStdAlnRep0/motif4
  Motif = Yale/Max:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562MaxStdAlnRep0/motif5

  Motif = Stanford/Mxi1_(AF4185):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mxi1af4185IggrabAlnRep0/motif1
  Motif = Stanford/Mxi1_(AF4185):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mxi1af4185IggrabAlnRep0/motif2
  Motif = Stanford/Mxi1_(AF4185):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mxi1af4185IggrabAlnRep0/motif3
  Motif = Stanford/Mxi1_(AF4185):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mxi1af4185IggrabAlnRep0/motif4
  Motif = Stanford/Mxi1_(AF4185):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Mxi1af4185IggrabAlnRep0/motif5

  Motif = Harvard/NELFe:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NelfeStdAlnRep0/motif1
  Motif = Harvard/NELFe:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NelfeStdAlnRep0/motif2
  Motif = Harvard/NELFe:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NelfeStdAlnRep0/motif3
  Motif = Harvard/NELFe:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NelfeStdAlnRep0/motif4
  Motif = Harvard/NELFe:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NelfeStdAlnRep0/motif5

  Motif = Yale/NF-E2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nfe2StdAlnRep0/motif1
  Motif = Yale/NF-E2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nfe2StdAlnRep0/motif2
  Motif = Yale/NF-E2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nfe2StdAlnRep0/motif3
  Motif = Yale/NF-E2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nfe2StdAlnRep0/motif4
  Motif = Yale/NF-E2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nfe2StdAlnRep0/motif5

  Motif = Stanford/NF-YA:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfyaStdAlnRep0/motif1
  Motif = Stanford/NF-YA:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfyaStdAlnRep0/motif2
  Motif = Stanford/NF-YA:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfyaStdAlnRep0/motif3
  Motif = Stanford/NF-YA:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfyaStdAlnRep0/motif4
  Motif = Stanford/NF-YA:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfyaStdAlnRep0/motif5

  Motif = Stanford/NF-YB:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfybStdAlnRep0/motif1
  Motif = Stanford/NF-YB:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfybStdAlnRep0/motif2
  Motif = Stanford/NF-YB:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfybStdAlnRep0/motif3
  Motif = Stanford/NF-YB:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfybStdAlnRep0/motif4
  Motif = Stanford/NF-YB:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562NfybStdAlnRep0/motif5

  Motif = Stanford/Nrf1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nrf1IggrabAlnRep0/motif1
  Motif = Stanford/Nrf1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nrf1IggrabAlnRep0/motif2
  Motif = Stanford/Nrf1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nrf1IggrabAlnRep0/motif3
  Motif = Stanford/Nrf1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nrf1IggrabAlnRep0/motif4
  Motif = Stanford/Nrf1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Nrf1IggrabAlnRep0/motif5

  Motif = Yale/Rad21:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rad21StdAlnRep0/motif1
  Motif = Yale/Rad21:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rad21StdAlnRep0/motif2
  Motif = Yale/Rad21:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rad21StdAlnRep0/motif3
  Motif = Yale/Rad21:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rad21StdAlnRep0/motif4
  Motif = Yale/Rad21:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rad21StdAlnRep0/motif5

  Motif = Stanford/RFX5_(200-401-194):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsHepg2Rfx5200401194IggrabAlnRep0/motif1
  Motif = Stanford/RFX5_(200-401-194):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsHepg2Rfx5200401194IggrabAlnRep0/motif2
  Motif = Stanford/RFX5_(200-401-194):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsHepg2Rfx5200401194IggrabAlnRep0/motif3
  Motif = Stanford/RFX5_(200-401-194):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsHepg2Rfx5200401194IggrabAlnRep0/motif4
  Motif = Stanford/RFX5_(200-401-194):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsHepg2Rfx5200401194IggrabAlnRep0/motif5

  Motif = Harvard/RPC155:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rpc155StdAlnRep0/motif1
  Motif = Harvard/RPC155:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rpc155StdAlnRep0/motif2
  Motif = Harvard/RPC155:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rpc155StdAlnRep0/motif3
  Motif = Harvard/RPC155:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rpc155StdAlnRep0/motif4
  Motif = Harvard/RPC155:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Rpc155StdAlnRep0/motif5

  Motif = Harvard/SIRT6:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Sirt6StdAlnRep0/motif1
  Motif = Harvard/SIRT6:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Sirt6StdAlnRep0/motif2
  Motif = Harvard/SIRT6:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Sirt6StdAlnRep0/motif3
  Motif = Harvard/SIRT6:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Sirt6StdAlnRep0/motif4
  Motif = Harvard/SIRT6:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Sirt6StdAlnRep0/motif5

  Motif = Stanford/SMC3_(ab9263):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Smc3ab9263IggrabAlnRep0/motif1
  Motif = Stanford/SMC3_(ab9263):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Smc3ab9263IggrabAlnRep0/motif2
  Motif = Stanford/SMC3_(ab9263):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Smc3ab9263IggrabAlnRep0/motif3
  Motif = Stanford/SMC3_(ab9263):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Smc3ab9263IggrabAlnRep0/motif4
  Motif = Stanford/SMC3_(ab9263):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Smc3ab9263IggrabAlnRep0/motif5

  Motif = Yale/STAT1/IFNa30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna30StdAlnRep0/motif1
  Motif = Yale/STAT1/IFNa30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna30StdAlnRep0/motif2
  Motif = Yale/STAT1/IFNa30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna30StdAlnRep0/motif3
  Motif = Yale/STAT1/IFNa30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna30StdAlnRep0/motif4
  Motif = Yale/STAT1/IFNa30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna30StdAlnRep0/motif5

  Motif = Yale/STAT1/IFNa6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna6hStdAlnRep0/motif1
  Motif = Yale/STAT1/IFNa6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna6hStdAlnRep0/motif2
  Motif = Yale/STAT1/IFNa6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna6hStdAlnRep0/motif3
  Motif = Yale/STAT1/IFNa6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna6hStdAlnRep0/motif4
  Motif = Yale/STAT1/IFNa6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifna6hStdAlnRep0/motif5

  Motif = Stanford/STAT1/IFNg30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng30StdAlnRep0/motif1
  Motif = Stanford/STAT1/IFNg30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng30StdAlnRep0/motif2
  Motif = Stanford/STAT1/IFNg30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng30StdAlnRep0/motif3
  Motif = Stanford/STAT1/IFNg30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng30StdAlnRep0/motif4
  Motif = Stanford/STAT1/IFNg30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng30StdAlnRep0/motif5

  Motif = Stanford/STAT1/IFNg6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng6hStdAlnRep0/motif1
  Motif = Stanford/STAT1/IFNg6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng6hStdAlnRep0/motif2
  Motif = Stanford/STAT1/IFNg6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng6hStdAlnRep0/motif3
  Motif = Stanford/STAT1/IFNg6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng6hStdAlnRep0/motif4
  Motif = Stanford/STAT1/IFNg6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat1Ifng6hStdAlnRep0/motif5

  Motif = Yale/STAT2/IFNa30:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna30StdAlnRep0/motif1
  Motif = Yale/STAT2/IFNa30:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna30StdAlnRep0/motif2
  Motif = Yale/STAT2/IFNa30:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna30StdAlnRep0/motif3
  Motif = Yale/STAT2/IFNa30:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna30StdAlnRep0/motif4
  Motif = Yale/STAT2/IFNa30:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna30StdAlnRep0/motif5

  Motif = Yale/STAT2/IFNa6h:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna6hStdAlnRep0/motif1
  Motif = Yale/STAT2/IFNa6h:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna6hStdAlnRep0/motif2
  Motif = Yale/STAT2/IFNa6h:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna6hStdAlnRep0/motif3
  Motif = Yale/STAT2/IFNa6h:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna6hStdAlnRep0/motif4
  Motif = Yale/STAT2/IFNa6h:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Stat2Ifna6hStdAlnRep0/motif5

  Motif = Stanford/TAL1_(SC-12984):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tal1sc12984IggmusAlnRep0/motif1
  Motif = Stanford/TAL1_(SC-12984):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tal1sc12984IggmusAlnRep0/motif2
  Motif = Stanford/TAL1_(SC-12984):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tal1sc12984IggmusAlnRep0/motif3
  Motif = Stanford/TAL1_(SC-12984):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tal1sc12984IggmusAlnRep0/motif4
  Motif = Stanford/TAL1_(SC-12984):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tal1sc12984IggmusAlnRep0/motif5

  Motif = Stanford/TBP:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562TbpIggmusAlnRep0/motif1
  Motif = Stanford/TBP:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562TbpIggmusAlnRep0/motif2
  Motif = Stanford/TBP:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562TbpIggmusAlnRep0/motif3
  Motif = Stanford/TBP:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562TbpIggmusAlnRep0/motif4
  Motif = Stanford/TBP:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562TbpIggmusAlnRep0/motif5

  Motif = Harvard/TFIIIC-110:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tf3c110StdAlnRep0/motif1
  Motif = Harvard/TFIIIC-110:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tf3c110StdAlnRep0/motif2
  Motif = Harvard/TFIIIC-110:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tf3c110StdAlnRep0/motif3
  Motif = Harvard/TFIIIC-110:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tf3c110StdAlnRep0/motif4
  Motif = Harvard/TFIIIC-110:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tf3c110StdAlnRep0/motif5

  Motif = USC/TR4:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tr4UcdAlnRep0/motif1
  Motif = USC/TR4:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tr4UcdAlnRep0/motif2
  Motif = USC/TR4:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tr4UcdAlnRep0/motif3
  Motif = USC/TR4:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tr4UcdAlnRep0/motif4
  Motif = USC/TR4:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Tr4UcdAlnRep0/motif5

  Motif = Stanford/USF2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Usf2IggrabAlnRep0/motif1
  Motif = Stanford/USF2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Usf2IggrabAlnRep0/motif2
  Motif = Stanford/USF2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Usf2IggrabAlnRep0/motif3
  Motif = Stanford/USF2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Usf2IggrabAlnRep0/motif4
  Motif = Stanford/USF2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Usf2IggrabAlnRep0/motif5

  Motif = USC/YY1:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Yy1UcdAlnRep0/motif1
  Motif = USC/YY1:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Yy1UcdAlnRep0/motif2
  Motif = USC/YY1:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Yy1UcdAlnRep0/motif3
  Motif = USC/YY1:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Yy1UcdAlnRep0/motif4
  Motif = USC/YY1:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Yy1UcdAlnRep0/motif5

  Motif = Stanford/Znf143_(16618-1-AP):1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsGm12878Znf143166181apStdAlnRep0/motif1
  Motif = Stanford/Znf143_(16618-1-AP):2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsGm12878Znf143166181apStdAlnRep0/motif2
  Motif = Stanford/Znf143_(16618-1-AP):3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsGm12878Znf143166181apStdAlnRep0/motif3
  Motif = Stanford/Znf143_(16618-1-AP):4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsGm12878Znf143166181apStdAlnRep0/motif4
  Motif = Stanford/Znf143_(16618-1-AP):5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsGm12878Znf143166181apStdAlnRep0/motif5

  Motif = USC/ZNF263:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Znf263UcdAlnRep0/motif1
  Motif = USC/ZNF263:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Znf263UcdAlnRep0/motif2
  Motif = USC/ZNF263:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Znf263UcdAlnRep0/motif3
  Motif = USC/ZNF263:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Znf263UcdAlnRep0/motif4
  Motif = USC/ZNF263:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeSydhTfbsK562Znf263UcdAlnRep0/motif5

  Motif = UChicago/eGFP-FOS:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EfosControlAlnRep0/motif1
  Motif = UChicago/eGFP-FOS:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EfosControlAlnRep0/motif2
  Motif = UChicago/eGFP-FOS:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EfosControlAlnRep0/motif3
  Motif = UChicago/eGFP-FOS:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EfosControlAlnRep0/motif4
  Motif = UChicago/eGFP-FOS:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EfosControlAlnRep0/motif5

  Motif = UChicago/eGFP-GATA2:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Egata2ControlAlnRep0/motif1
  Motif = UChicago/eGFP-GATA2:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Egata2ControlAlnRep0/motif2
  Motif = UChicago/eGFP-GATA2:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Egata2ControlAlnRep0/motif3
  Motif = UChicago/eGFP-GATA2:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Egata2ControlAlnRep0/motif4
  Motif = UChicago/eGFP-GATA2:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Egata2ControlAlnRep0/motif5

  Motif = UChicago/eGFP-HDAC8:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Ehdac8ControlAlnRep0/motif1
  Motif = UChicago/eGFP-HDAC8:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Ehdac8ControlAlnRep0/motif2
  Motif = UChicago/eGFP-HDAC8:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Ehdac8ControlAlnRep0/motif3
  Motif = UChicago/eGFP-HDAC8:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Ehdac8ControlAlnRep0/motif4
  Motif = UChicago/eGFP-HDAC8:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562Ehdac8ControlAlnRep0/motif5

  Motif = UChicago/eGFP-JunB:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjunbControlAlnRep0/motif1
  Motif = UChicago/eGFP-JunB:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjunbControlAlnRep0/motif2
  Motif = UChicago/eGFP-JunB:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjunbControlAlnRep0/motif3
  Motif = UChicago/eGFP-JunB:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjunbControlAlnRep0/motif4
  Motif = UChicago/eGFP-JunB:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjunbControlAlnRep0/motif5

  Motif = UChicago/eGFP-JunD:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjundControlAlnRep0/motif1
  Motif = UChicago/eGFP-JunD:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjundControlAlnRep0/motif2
  Motif = UChicago/eGFP-JunD:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjundControlAlnRep0/motif3
  Motif = UChicago/eGFP-JunD:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjundControlAlnRep0/motif4
  Motif = UChicago/eGFP-JunD:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUchicagoTfbsK562EjundControlAlnRep0/motif5

  Motif = UW/CTCF:1 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif1
  Motif = UW/CTCF:2 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif2
  Motif = UW/CTCF:3 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif3
  Motif = UW/CTCF:4 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif4
  Motif = UW/CTCF:5 K562_ChIP-seq_motifs/spp.optimal.wgEncodeUwTfbsK562CtcfStdAlnRep0/motif5

  Sensitivity = 0.8
</Motifs>

#
# Scope of the analysis, limiting the set of motif complexes or datasets considered
#

<Scope>
  Dataset = Broad/CTCF
  Motif1 = Broad/CTCF:1 Broad/CTCF:2 Broad/CTCF:3 Broad/CTCF:4 Broad/CTCF:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/ATF3
  Motif1 = HudsonAlpha/ATF3:1 HudsonAlpha/ATF3:2 HudsonAlpha/ATF3:3 HudsonAlpha/ATF3:4 HudsonAlpha/ATF3:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/BCLAF1_(SC-101388)
  Motif1 = HudsonAlpha/BCLAF1_(SC-101388):1 HudsonAlpha/BCLAF1_(SC-101388):2 HudsonAlpha/BCLAF1_(SC-101388):3 HudsonAlpha/BCLAF1_(SC-101388):4 HudsonAlpha/BCLAF1_(SC-101388):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/CTCF_(SC-5916)
  Motif1 = HudsonAlpha/CTCF_(SC-5916):1 HudsonAlpha/CTCF_(SC-5916):2 HudsonAlpha/CTCF_(SC-5916):3 HudsonAlpha/CTCF_(SC-5916):4 HudsonAlpha/CTCF_(SC-5916):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/CTCFL_(SC-98982)
  Motif1 = HudsonAlpha/CTCFL_(SC-98982):1 HudsonAlpha/CTCFL_(SC-98982):2 HudsonAlpha/CTCFL_(SC-98982):3 HudsonAlpha/CTCFL_(SC-98982):4 HudsonAlpha/CTCFL_(SC-98982):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/E2F6
  Motif1 = HudsonAlpha/E2F6:1 HudsonAlpha/E2F6:2 HudsonAlpha/E2F6:3 HudsonAlpha/E2F6:4 HudsonAlpha/E2F6:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/Egr-1
  Motif1 = HudsonAlpha/Egr-1:1 HudsonAlpha/Egr-1:2 HudsonAlpha/Egr-1:3 HudsonAlpha/Egr-1:4 HudsonAlpha/Egr-1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/ELF1_(SC-631)
  Motif1 = HudsonAlpha/ELF1_(SC-631):1 HudsonAlpha/ELF1_(SC-631):2 HudsonAlpha/ELF1_(SC-631):3 HudsonAlpha/ELF1_(SC-631):4 HudsonAlpha/ELF1_(SC-631):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/ETS1
  Motif1 = HudsonAlpha/ETS1:1 HudsonAlpha/ETS1:2 HudsonAlpha/ETS1:3 HudsonAlpha/ETS1:4 HudsonAlpha/ETS1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/FOSL1_(SC-183)
  Motif1 = HudsonAlpha/FOSL1_(SC-183):1 HudsonAlpha/FOSL1_(SC-183):2 HudsonAlpha/FOSL1_(SC-183):3 HudsonAlpha/FOSL1_(SC-183):4 HudsonAlpha/FOSL1_(SC-183):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/GABP
  Motif1 = HudsonAlpha/GABP:1 HudsonAlpha/GABP:2 HudsonAlpha/GABP:3 HudsonAlpha/GABP:4 HudsonAlpha/GABP:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/GATA2_(SC-267)
  Motif1 = HudsonAlpha/GATA2_(SC-267):1 HudsonAlpha/GATA2_(SC-267):2 HudsonAlpha/GATA2_(SC-267):3 HudsonAlpha/GATA2_(SC-267):4 HudsonAlpha/GATA2_(SC-267):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/HDAC2_(SC-6296)
  Motif1 = HudsonAlpha/HDAC2_(SC-6296):1 HudsonAlpha/HDAC2_(SC-6296):2 HudsonAlpha/HDAC2_(SC-6296):3 HudsonAlpha/HDAC2_(SC-6296):4 HudsonAlpha/HDAC2_(SC-6296):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/Max
  Motif1 = HudsonAlpha/Max:1 HudsonAlpha/Max:2 HudsonAlpha/Max:3 HudsonAlpha/Max:4 HudsonAlpha/Max:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/MEF2A
  Motif1 = HudsonAlpha/MEF2A:1 HudsonAlpha/MEF2A:2 HudsonAlpha/MEF2A:3 HudsonAlpha/MEF2A:4 HudsonAlpha/MEF2A:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/NRSF
  Motif1 = HudsonAlpha/NRSF:1 HudsonAlpha/NRSF:2 HudsonAlpha/NRSF:3 HudsonAlpha/NRSF:4 HudsonAlpha/NRSF:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/PU.1
  Motif1 = HudsonAlpha/PU.1:1 HudsonAlpha/PU.1:2 HudsonAlpha/PU.1:3 HudsonAlpha/PU.1:4 HudsonAlpha/PU.1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/Rad21
  Motif1 = HudsonAlpha/Rad21:1 HudsonAlpha/Rad21:2 HudsonAlpha/Rad21:3 HudsonAlpha/Rad21:4 HudsonAlpha/Rad21:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/Sin3Ak-20
  Motif1 = HudsonAlpha/Sin3Ak-20:1 HudsonAlpha/Sin3Ak-20:2 HudsonAlpha/Sin3Ak-20:3 HudsonAlpha/Sin3Ak-20:4 HudsonAlpha/Sin3Ak-20:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/SIX5
  Motif1 = HudsonAlpha/SIX5:1 HudsonAlpha/SIX5:2 HudsonAlpha/SIX5:3 HudsonAlpha/SIX5:4 HudsonAlpha/SIX5:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/SP1
  Motif1 = HudsonAlpha/SP1:1 HudsonAlpha/SP1:2 HudsonAlpha/SP1:3 HudsonAlpha/SP1:4 HudsonAlpha/SP1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/SP2_(SC-643)
  Motif1 = HudsonAlpha/SP2_(SC-643):1 HudsonAlpha/SP2_(SC-643):2 HudsonAlpha/SP2_(SC-643):3 HudsonAlpha/SP2_(SC-643):4 HudsonAlpha/SP2_(SC-643):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/SRF
  Motif1 = HudsonAlpha/SRF:1 HudsonAlpha/SRF:2 HudsonAlpha/SRF:3 HudsonAlpha/SRF:4 HudsonAlpha/SRF:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/TAF1
  Motif1 = HudsonAlpha/TAF1:1 HudsonAlpha/TAF1:2 HudsonAlpha/TAF1:3 HudsonAlpha/TAF1:4 HudsonAlpha/TAF1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/TAF7_(SC-101167)
  Motif1 = HudsonAlpha/TAF7_(SC-101167):1 HudsonAlpha/TAF7_(SC-101167):2 HudsonAlpha/TAF7_(SC-101167):3 HudsonAlpha/TAF7_(SC-101167):4 HudsonAlpha/TAF7_(SC-101167):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/THAP1_(SC-98174)
  Motif1 = HudsonAlpha/THAP1_(SC-98174):1 HudsonAlpha/THAP1_(SC-98174):2 HudsonAlpha/THAP1_(SC-98174):3 HudsonAlpha/THAP1_(SC-98174):4 HudsonAlpha/THAP1_(SC-98174):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/USF-1
  Motif1 = HudsonAlpha/USF-1:1 HudsonAlpha/USF-1:2 HudsonAlpha/USF-1:3 HudsonAlpha/USF-1:4 HudsonAlpha/USF-1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/YY1
  Motif1 = HudsonAlpha/YY1:1 HudsonAlpha/YY1:2 HudsonAlpha/YY1:3 HudsonAlpha/YY1:4 HudsonAlpha/YY1:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/YY1_(SC-281)
  Motif1 = HudsonAlpha/YY1_(SC-281):1 HudsonAlpha/YY1_(SC-281):2 HudsonAlpha/YY1_(SC-281):3 HudsonAlpha/YY1_(SC-281):4 HudsonAlpha/YY1_(SC-281):5
</Scope>

<Scope>
  Dataset = HudsonAlpha/ZBTB33
  Motif1 = HudsonAlpha/ZBTB33:1 HudsonAlpha/ZBTB33:2 HudsonAlpha/ZBTB33:3 HudsonAlpha/ZBTB33:4 HudsonAlpha/ZBTB33:5
</Scope>

<Scope>
  Dataset = HudsonAlpha/ZBTB7A_(SC-34508)
  Motif1 = HudsonAlpha/ZBTB7A_(SC-34508):1 HudsonAlpha/ZBTB7A_(SC-34508):2 HudsonAlpha/ZBTB7A_(SC-34508):3 HudsonAlpha/ZBTB7A_(SC-34508):4 HudsonAlpha/ZBTB7A_(SC-34508):5
</Scope>

<Scope>
  Dataset = Harvard/ATF3
  Motif1 = Harvard/ATF3:1 Harvard/ATF3:2 Harvard/ATF3:3 Harvard/ATF3:4 Harvard/ATF3:5
</Scope>

<Scope>
  Dataset = Harvard/BDP1
  Motif1 = Harvard/BDP1:1 Harvard/BDP1:2 Harvard/BDP1:3 Harvard/BDP1:4 Harvard/BDP1:5
</Scope>

<Scope>
  Dataset = Stanford/BHLHE40_(NB100-1800)
  Motif1 = Stanford/BHLHE40_(NB100-1800):1 Stanford/BHLHE40_(NB100-1800):2 Stanford/BHLHE40_(NB100-1800):3 Stanford/BHLHE40_(NB100-1800):4 Stanford/BHLHE40_(NB100-1800):5
</Scope>

<Scope>
  Dataset = Harvard/BRF1
  Motif1 = Harvard/BRF1:1 Harvard/BRF1:2 Harvard/BRF1:3 Harvard/BRF1:4 Harvard/BRF1:5
</Scope>

<Scope>
  Dataset = Harvard/BRF2
  Motif1 = Harvard/BRF2:1 Harvard/BRF2:2 Harvard/BRF2:3 Harvard/BRF2:4 Harvard/BRF2:5
</Scope>

<Scope>
  Dataset = Stanford/Brg1
  Motif1 = Stanford/Brg1:1 Stanford/Brg1:2 Stanford/Brg1:3 Stanford/Brg1:4 Stanford/Brg1:5
</Scope>

<Scope>
  Dataset = Harvard/CCNT2
  Motif1 = Harvard/CCNT2:1 Harvard/CCNT2:2 Harvard/CCNT2:3 Harvard/CCNT2:4 Harvard/CCNT2:5
</Scope>

<Scope>
  Dataset = Stanford/CEBPB
  Motif1 = Stanford/CEBPB:1 Stanford/CEBPB:2 Stanford/CEBPB:3 Stanford/CEBPB:4 Stanford/CEBPB:5
</Scope>

<Scope>
  Dataset = Yale/c-Fos
  Motif1 = Yale/c-Fos:1 Yale/c-Fos:2 Yale/c-Fos:3 Yale/c-Fos:4 Yale/c-Fos:5
</Scope>

<Scope>
  Dataset = Stanford/CHD2_(AB68301)
  Motif1 = Stanford/CHD2_(AB68301):1 Stanford/CHD2_(AB68301):2 Stanford/CHD2_(AB68301):3 Stanford/CHD2_(AB68301):4 Stanford/CHD2_(AB68301):5
</Scope>

<Scope>
  Dataset = Yale/c-Jun/IFNa6h
  Motif1 = Yale/c-Jun/IFNa6h:1 Yale/c-Jun/IFNa6h:2 Yale/c-Jun/IFNa6h:3 Yale/c-Jun/IFNa6h:4 Yale/c-Jun/IFNa6h:5
</Scope>

<Scope>
  Dataset = Yale/c-Jun/IFNg30
  Motif1 = Yale/c-Jun/IFNg30:1 Yale/c-Jun/IFNg30:2 Yale/c-Jun/IFNg30:3 Yale/c-Jun/IFNg30:4 Yale/c-Jun/IFNg30:5
</Scope>

<Scope>
  Dataset = Yale/c-Jun/IFNg6h
  Motif1 = Yale/c-Jun/IFNg6h:1 Yale/c-Jun/IFNg6h:2 Yale/c-Jun/IFNg6h:3 Yale/c-Jun/IFNg6h:4 Yale/c-Jun/IFNg6h:5
</Scope>

<Scope>
  Dataset = Yale/c-Jun
  Motif1 = Yale/c-Jun:1 Yale/c-Jun:2 Yale/c-Jun:3 Yale/c-Jun:4 Yale/c-Jun:5
</Scope>

<Scope>
  Dataset = Yale/c-Myc/IFNa30
  Motif1 = Yale/c-Myc/IFNa30:1 Yale/c-Myc/IFNa30:2 Yale/c-Myc/IFNa30:3 Yale/c-Myc/IFNa30:4 Yale/c-Myc/IFNa30:5
</Scope>

<Scope>
  Dataset = Yale/c-Myc/IFNa6h
  Motif1 = Yale/c-Myc/IFNa6h:1 Yale/c-Myc/IFNa6h:2 Yale/c-Myc/IFNa6h:3 Yale/c-Myc/IFNa6h:4 Yale/c-Myc/IFNa6h:5
</Scope>

<Scope>
  Dataset = Stanford/c-Myc/IFNg30
  Motif1 = Stanford/c-Myc/IFNg30:1 Stanford/c-Myc/IFNg30:2 Stanford/c-Myc/IFNg30:3 Stanford/c-Myc/IFNg30:4 Stanford/c-Myc/IFNg30:5
</Scope>

<Scope>
  Dataset = Yale/c-Myc/IFNg6h
  Motif1 = Yale/c-Myc/IFNg6h:1 Yale/c-Myc/IFNg6h:2 Yale/c-Myc/IFNg6h:3 Yale/c-Myc/IFNg6h:4 Yale/c-Myc/IFNg6h:5
</Scope>

<Scope>
  Dataset = Yale/c-Myc
  Motif1 = Yale/c-Myc:1 Yale/c-Myc:2 Yale/c-Myc:3 Yale/c-Myc:4 Yale/c-Myc:5
</Scope>

<Scope>
  Dataset = USC/E2F4
  Motif1 = USC/E2F4:1 USC/E2F4:2 USC/E2F4:3 USC/E2F4:4 USC/E2F4:5
</Scope>

<Scope>
  Dataset = USC/E2F6
  Motif1 = USC/E2F6:1 USC/E2F6:2 USC/E2F6:3 USC/E2F6:4 USC/E2F6:5
</Scope>

<Scope>
  Dataset = USC/GATA-1
  Motif1 = USC/GATA-1:1 USC/GATA-1:2 USC/GATA-1:3 USC/GATA-1:4 USC/GATA-1:5
</Scope>

<Scope>
  Dataset = USC/GATA-2
  Motif1 = USC/GATA-2:1 USC/GATA-2:2 USC/GATA-2:3 USC/GATA-2:4 USC/GATA-2:5
</Scope>

<Scope>
  Dataset = Harvard/GTF2B
  Motif1 = Harvard/GTF2B:1 Harvard/GTF2B:2 Harvard/GTF2B:3 Harvard/GTF2B:4 Harvard/GTF2B:5
</Scope>

<Scope>
  Dataset = Stanford/GTF2F1_(AB28179)
  Motif1 = Stanford/GTF2F1_(AB28179):1 Stanford/GTF2F1_(AB28179):2 Stanford/GTF2F1_(AB28179):3 Stanford/GTF2F1_(AB28179):4 Stanford/GTF2F1_(AB28179):5
</Scope>

<Scope>
  Dataset = Harvard/HMGN3
  Motif1 = Harvard/HMGN3:1 Harvard/HMGN3:2 Harvard/HMGN3:3 Harvard/HMGN3:4 Harvard/HMGN3:5
</Scope>

<Scope>
  Dataset = Stanford/Ini1
  Motif1 = Stanford/Ini1:1 Stanford/Ini1:2 Stanford/Ini1:3 Stanford/Ini1:4 Stanford/Ini1:5
</Scope>

<Scope>
  Dataset = Stanford/IRF1/IFNa30
  Motif1 = Stanford/IRF1/IFNa30:1 Stanford/IRF1/IFNa30:2 Stanford/IRF1/IFNa30:3 Stanford/IRF1/IFNa30:4 Stanford/IRF1/IFNa30:5
</Scope>

<Scope>
  Dataset = Stanford/IRF1/IFNg6h
  Motif1 = Stanford/IRF1/IFNg6h:1 Stanford/IRF1/IFNg6h:2 Stanford/IRF1/IFNg6h:3 Stanford/IRF1/IFNg6h:4 Stanford/IRF1/IFNg6h:5
</Scope>

<Scope>
  Dataset = Stanford/JunD
  Motif1 = Stanford/JunD:1 Stanford/JunD:2 Stanford/JunD:3 Stanford/JunD:4 Stanford/JunD:5
</Scope>

<Scope>
  Dataset = Stanford/MafK_(ab50322)
  Motif1 = Stanford/MafK_(ab50322):1 Stanford/MafK_(ab50322):2 Stanford/MafK_(ab50322):3 Stanford/MafK_(ab50322):4 Stanford/MafK_(ab50322):5
</Scope>

<Scope>
  Dataset = Yale/Max
  Motif1 = Yale/Max:1 Yale/Max:2 Yale/Max:3 Yale/Max:4 Yale/Max:5
</Scope>

<Scope>
  Dataset = Stanford/Mxi1_(AF4185)
  Motif1 = Stanford/Mxi1_(AF4185):1 Stanford/Mxi1_(AF4185):2 Stanford/Mxi1_(AF4185):3 Stanford/Mxi1_(AF4185):4 Stanford/Mxi1_(AF4185):5
</Scope>

<Scope>
  Dataset = Harvard/NELFe
  Motif1 = Harvard/NELFe:1 Harvard/NELFe:2 Harvard/NELFe:3 Harvard/NELFe:4 Harvard/NELFe:5
</Scope>

<Scope>
  Dataset = Yale/NF-E2
  Motif1 = Yale/NF-E2:1 Yale/NF-E2:2 Yale/NF-E2:3 Yale/NF-E2:4 Yale/NF-E2:5
</Scope>

<Scope>
  Dataset = Stanford/NF-YA
  Motif1 = Stanford/NF-YA:1 Stanford/NF-YA:2 Stanford/NF-YA:3 Stanford/NF-YA:4 Stanford/NF-YA:5
</Scope>

<Scope>
  Dataset = Stanford/NF-YB
  Motif1 = Stanford/NF-YB:1 Stanford/NF-YB:2 Stanford/NF-YB:3 Stanford/NF-YB:4 Stanford/NF-YB:5
</Scope>

<Scope>
  Dataset = Stanford/Nrf1
  Motif1 = Stanford/Nrf1:1 Stanford/Nrf1:2 Stanford/Nrf1:3 Stanford/Nrf1:4 Stanford/Nrf1:5
</Scope>

<Scope>
  Dataset = Yale/Rad21
  Motif1 = Yale/Rad21:1 Yale/Rad21:2 Yale/Rad21:3 Yale/Rad21:4 Yale/Rad21:5
</Scope>

<Scope>
  Dataset = Stanford/RFX5_(200-401-194)
  Motif1 = Stanford/RFX5_(200-401-194):1 Stanford/RFX5_(200-401-194):2 Stanford/RFX5_(200-401-194):3 Stanford/RFX5_(200-401-194):4 Stanford/RFX5_(200-401-194):5
</Scope>

<Scope>
  Dataset = Harvard/RPC155
  Motif1 = Harvard/RPC155:1 Harvard/RPC155:2 Harvard/RPC155:3 Harvard/RPC155:4 Harvard/RPC155:5
</Scope>

<Scope>
  Dataset = Harvard/SIRT6
  Motif1 = Harvard/SIRT6:1 Harvard/SIRT6:2 Harvard/SIRT6:3 Harvard/SIRT6:4 Harvard/SIRT6:5
</Scope>

<Scope>
  Dataset = Stanford/SMC3_(ab9263)
  Motif1 = Stanford/SMC3_(ab9263):1 Stanford/SMC3_(ab9263):2 Stanford/SMC3_(ab9263):3 Stanford/SMC3_(ab9263):4 Stanford/SMC3_(ab9263):5
</Scope>

<Scope>
  Dataset = Yale/STAT1/IFNa30
  Motif1 = Yale/STAT1/IFNa30:1 Yale/STAT1/IFNa30:2 Yale/STAT1/IFNa30:3 Yale/STAT1/IFNa30:4 Yale/STAT1/IFNa30:5
</Scope>

<Scope>
  Dataset = Yale/STAT1/IFNa6h
  Motif1 = Yale/STAT1/IFNa6h:1 Yale/STAT1/IFNa6h:2 Yale/STAT1/IFNa6h:3 Yale/STAT1/IFNa6h:4 Yale/STAT1/IFNa6h:5
</Scope>

<Scope>
  Dataset = Stanford/STAT1/IFNg30
  Motif1 = Stanford/STAT1/IFNg30:1 Stanford/STAT1/IFNg30:2 Stanford/STAT1/IFNg30:3 Stanford/STAT1/IFNg30:4 Stanford/STAT1/IFNg30:5
</Scope>

<Scope>
  Dataset = Stanford/STAT1/IFNg6h
  Motif1 = Stanford/STAT1/IFNg6h:1 Stanford/STAT1/IFNg6h:2 Stanford/STAT1/IFNg6h:3 Stanford/STAT1/IFNg6h:4 Stanford/STAT1/IFNg6h:5
</Scope>

<Scope>
  Dataset = Yale/STAT2/IFNa30
  Motif1 = Yale/STAT2/IFNa30:1 Yale/STAT2/IFNa30:2 Yale/STAT2/IFNa30:3 Yale/STAT2/IFNa30:4 Yale/STAT2/IFNa30:5
</Scope>

<Scope>
  Dataset = Yale/STAT2/IFNa6h
  Motif1 = Yale/STAT2/IFNa6h:1 Yale/STAT2/IFNa6h:2 Yale/STAT2/IFNa6h:3 Yale/STAT2/IFNa6h:4 Yale/STAT2/IFNa6h:5
</Scope>

<Scope>
  Dataset = Stanford/TAL1_(SC-12984)
  Motif1 = Stanford/TAL1_(SC-12984):1 Stanford/TAL1_(SC-12984):2 Stanford/TAL1_(SC-12984):3 Stanford/TAL1_(SC-12984):4 Stanford/TAL1_(SC-12984):5
</Scope>

<Scope>
  Dataset = Stanford/TBP
  Motif1 = Stanford/TBP:1 Stanford/TBP:2 Stanford/TBP:3 Stanford/TBP:4 Stanford/TBP:5
</Scope>

<Scope>
  Dataset = Harvard/TFIIIC-110
  Motif1 = Harvard/TFIIIC-110:1 Harvard/TFIIIC-110:2 Harvard/TFIIIC-110:3 Harvard/TFIIIC-110:4 Harvard/TFIIIC-110:5
</Scope>

<Scope>
  Dataset = USC/TR4
  Motif1 = USC/TR4:1 USC/TR4:2 USC/TR4:3 USC/TR4:4 USC/TR4:5
</Scope>

<Scope>
  Dataset = Stanford/USF2
  Motif1 = Stanford/USF2:1 Stanford/USF2:2 Stanford/USF2:3 Stanford/USF2:4 Stanford/USF2:5
</Scope>

<Scope>
  Dataset = USC/YY1
  Motif1 = USC/YY1:1 USC/YY1:2 USC/YY1:3 USC/YY1:4 USC/YY1:5
</Scope>

<Scope>
  Dataset = Stanford/Znf143_(16618-1-AP)
  Motif1 = Stanford/Znf143_(16618-1-AP):1 Stanford/Znf143_(16618-1-AP):2 Stanford/Znf143_(16618-1-AP):3 Stanford/Znf143_(16618-1-AP):4 Stanford/Znf143_(16618-1-AP):5
</Scope>

<Scope>
  Dataset = USC/ZNF263
  Motif1 = USC/ZNF263:1 USC/ZNF263:2 USC/ZNF263:3 USC/ZNF263:4 USC/ZNF263:5
</Scope>

<Scope>
  Dataset = UChicago/eGFP-FOS
  Motif1 = UChicago/eGFP-FOS:1 UChicago/eGFP-FOS:2 UChicago/eGFP-FOS:3 UChicago/eGFP-FOS:4 UChicago/eGFP-FOS:5
</Scope>

<Scope>
  Dataset = UChicago/eGFP-GATA2
  Motif1 = UChicago/eGFP-GATA2:1 UChicago/eGFP-GATA2:2 UChicago/eGFP-GATA2:3 UChicago/eGFP-GATA2:4 UChicago/eGFP-GATA2:5
</Scope>

<Scope>
  Dataset = UChicago/eGFP-HDAC8
  Motif1 = UChicago/eGFP-HDAC8:1 UChicago/eGFP-HDAC8:2 UChicago/eGFP-HDAC8:3 UChicago/eGFP-HDAC8:4 UChicago/eGFP-HDAC8:5
</Scope>

<Scope>
  Dataset = UChicago/eGFP-JunB
  Motif1 = UChicago/eGFP-JunB:1 UChicago/eGFP-JunB:2 UChicago/eGFP-JunB:3 UChicago/eGFP-JunB:4 UChicago/eGFP-JunB:5
</Scope>

<Scope>
  Dataset = UChicago/eGFP-JunD
  Motif1 = UChicago/eGFP-JunD:1 UChicago/eGFP-JunD:2 UChicago/eGFP-JunD:3 UChicago/eGFP-JunD:4 UChicago/eGFP-JunD:5
</Scope>

<Scope>
  Dataset = UW/CTCF
  Motif1 = UW/CTCF:1 UW/CTCF:2 UW/CTCF:3 UW/CTCF:4 UW/CTCF:5
</Scope>

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
  FoldChangeThreshold = 2.0
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
