import FWCore.ParameterSet.Config as cms

process = cms.Process("TTBSM")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( 
  wantSummary = cms.untracked.bool(True), 
  allowUnscheduled = cms.untracked.bool(True)
  )

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

##########################################################################################
useMiniAOD = True

# AOD
pfcandidates          = 'particleFlow'
chsstring             = 'pfNoPileUpJME'
genjetparticles       = 'genParticles'
importantgenparticles = 'genParticles'
tracks                = 'generalTracks'
vertices              = 'offlinePrimaryVertices'
mergedvertices        = 'inclusiveMergedVertices' 
mergedvertices2       = '' 
primaryvertices       = 'offlinePrimaryVertices'

#miniAOD
if useMiniAOD:
  pfcandidates          = 'packedPFCandidates'
  genjetparticles       = 'packedGenParticles'
  importantgenparticles = 'prunedGenParticles'
  tracks                = 'unpackedTracksAndVertices'
  vertices              = 'unpackedTracksAndVertices'
  mergedvertices        = 'unpackedTracksAndVertices'
  mergedvertices2       = 'secondary'
  primaryvertices       = 'offlineSlimmedPrimaryVertices'


print 'useMiniAOD : '+str(useMiniAOD)
print 'pfcandidates          : '+pfcandidates         
print 'genjetparticles       : '+genjetparticles      
print 'importantgenparticles : '+importantgenparticles
print 'tracks                : '+tracks               
print 'vertices              : '+vertices             
print 'mergedvertices        : '+mergedvertices       
print 'mergedvertices2       : '+mergedvertices2      
print 'primaryvertices       : '+primaryvertices 
##########################################################################################
# INPUT
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring(
            'root://cmsxrootd-site.fnal.gov//store/user/jdolen/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/140707_143029/0000/miniAOD-prod_PAT_1.root'
      	)
                                )
##########################################################################################
# OUTPUT
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("ttbsm_ultraslim.root"),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_ttbsmAna*_*_*'
                                                                      #, 'keep *_goodPatJetsCA8PrunedPF_*_*'
                                                                      #, 'keep *_goodPatJetsCATopTagPF_*_*'
                                                                      #, 'keep recoPFJets_*_*_*'
                                                                      ) 
                               )
##########################################################################################
# SETUP
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START70_V6::All'

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.fixedGridRhoFastjetAll.pfCandidatesTag = pfcandidates


from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.caTopTaggers_cff import *

##########################################################################################

## Geometry and Detector Conditions (needed for a few patTuple production steps)
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_R_42_V12::All'
#process.load("Configuration.StandardSequences.MagneticField_cff")

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
from Analysis.BoostedTopAnalysis.CATopTagParams_cfi import caTopTagParams
from Analysis.BoostedTopAnalysis.BoostedTopWTagParams_cfi import boostedTopWTagParams

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring('HLT_HT750_v*'),           # provide list of HLT paths (or patterns) you want
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True)    # throw exception on unknown path names
)

##########################################################################################
# GENJETS
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag(genjetparticles))

process.ca8GenJetsPruned = ca4GenJets.clone(
  rParam = cms.double(0.8),
  src = cms.InputTag(genjetparticles),
  usePruning = cms.bool(True),
  zcut = cms.double(0.1),
  rcut_factor = cms.double(0.5),
  nFilt = cms.int32(2),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )

##########################################################################################
# RECO

from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment


if useMiniAOD:
  process.chs = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("fromPV"))
  chsstring = 'chs'
process.ak4PFJets.src = pfcandidates
process.ak5PFJets.src = pfcandidates
process.ca4PFJets.src = pfcandidates
process.kt6PFJets.src = pfcandidates

# kt6PFJets = kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True, src = cms.InputTag(pfcandidates) ) 

process.ca8PFJets  = process.ca4PFJets.clone(rParam = 0.8,  doAreaFastjet = True)
process.ak15PFJets = process.ak4PFJets.clone(rParam = 1.5,  doAreaFastjet = True)

process.ca8CHSJets  = process.ca8PFJets.clone (src = chsstring)
process.ak15CHSJets = process.ak15PFJets.clone(src = chsstring)

process.load('RecoJets.JetProducers.ak4PFJetsPruned_cfi')
process.ca8PFJetsPruned  = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), doAreaFastjet = True, src = pfcandidates)
process.ca8CHSJetsPruned = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), doAreaFastjet = True, src = chsstring)

# CATopJet PF Jets with adjacency 
#process.cmsTopTagCHS = cmsTopTagPFJetsCHS.clone()
process.cmsTopTagCHS = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag(chsstring),
                           doAreaFastjet = cms.bool(True),
                           doRhoFastjet = cms.bool(False),
                           jetPtMin = cms.double(100.0)
                           ),
    AnomalousCellParameters,
    CATopJetParameters.clone( jetCollInstanceName = cms.string("SubJets"), 
                              verbose = cms.bool(False),
                              algorithm = cms.int32(1), # 0 = KT, 1 = CA, 2 = anti-KT
                              tagAlgo = cms.int32(0), #0=legacy top
                              useAdjacency = cms.int32(2), # modified adjacency
                              centralEtaCut = cms.double(2.5), # eta for defining "central" jets
                              sumEtBins = cms.vdouble(0,1600,2600), # sumEt bins over which cuts vary. vector={bin 0 lower bound, bin 1 lower bound, ...}
                              rBins = cms.vdouble(0.8,0.8,0.8), # Jet distance paramter R. R values depend on sumEt bins.
                              ptFracBins = cms.vdouble(0.05,0.05,0.05), # minimum fraction of central jet pt for subjets (deltap)
                              deltarBins = cms.vdouble(0.19,0.19,0.19), # Applicable only if useAdjacency=1. deltar adjacency values for each sumEtBin
                              nCellBins = cms.vdouble(1.9,1.9,1.9), 
                            ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )



process.CATopTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("cmsTopTagCHS"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )


###############################################
# PAT JETS
from PhysicsTools.PatAlgos.tools.jetTools import *

process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

# fix JTA (see https://github.com/cms-sw/cmssw/tree/CMSSW_7_0_X/RecoJets/JetAssociationProducers/python)

process.ak5JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJets")
process.ak5JetTracksAssociatorAtVertexPF.tracks = cms.InputTag(tracks)
process.impactParameterTagInfos.primaryVertex = cms.InputTag(vertices)
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag(mergedvertices,mergedvertices2,"")
process.combinedSecondaryVertex.trackMultiplicityMin = 1


# patJetsCA8PF
addJetCollection(
    process,
    labelName = 'CA8PF',
    jetSource = cms.InputTag('ca8PFJets'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8PF.matched=importantgenparticles
process.patJetCorrFactorsCA8PF.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8PF.matched = 'ca8GenJets'#'slimmedGenJets'
process.patJetPartons.particles = importantgenparticles
process.jetTracksAssociatorAtVertexCA8PF=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8PFJets'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8PF.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8PF.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8PF= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8PF.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8PF.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8PF.jetTagComputer = cms.string('combinedSecondaryVertexCA8PF')

#patJetsCA8CHS
addJetCollection(
    process,
    labelName = 'CA8CHS',
    jetSource = cms.InputTag('ca8CHSJets'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8CHS.matched=importantgenparticles
process.patJetCorrFactorsCA8CHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHS.matched = 'ca8GenJets'#'slimmedGenJets'
process.jetTracksAssociatorAtVertexCA8CHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8CHSJets'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8CHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8CHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8CHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8CHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8CHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8CHS.jetTagComputer = cms.string('combinedSecondaryVertexCA8CHS')


# patJetsCA8CHSpruned
addJetCollection(
    process,
    labelName = 'CA8CHSpruned',
    jetSource = cms.InputTag('ca8CHSJetsPruned'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8CHSpruned.matched=importantgenparticles
process.patJetCorrFactorsCA8CHSpruned.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHSpruned.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCA8CHSpruned.matched = importantgenparticles
process.jetTracksAssociatorAtVertexCA8CHSpruned=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8CHSJetsPruned'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8CHSpruned.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8CHSpruned.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8CHSpruned= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8CHSpruned.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8CHSpruned.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8CHSpruned.jetTagComputer = cms.string('combinedSecondaryVertexCA8CHSpruned')

addJetCollection(
    process,
    labelName = 'CA8CHSprunedSubjets',
    jetSource = cms.InputTag('ca8CHSJetsPruned','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )

process.patJetPartonMatchCA8CHSprunedSubjets.matched=importantgenparticles
process.patJetCorrFactorsCA8CHSprunedSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHSprunedSubjets.matched = 'ca8GenJetsPruned'#slimmedGenJets'
process.patJetPartonMatchCA8CHSprunedSubjets.matched = importantgenparticles

process.patJetsCA8CHSprunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCA8CHSpruned" ),
    subjetSrc=cms.InputTag("patJetsCA8CHSprunedSubjets")
      )


# patJetsCMSTopTagCHS
addJetCollection(
    process,
    labelName = 'CMSTopTagCHS',
    jetSource = cms.InputTag('cmsTopTagCHS'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
        # btagInfos = [
        # 'CATopTagInfos'
        #  ]   
    )
process.patJetPartonMatchCMSTopTagCHS.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagCHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagCHS.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHS.matched = importantgenparticles
process.jetTracksAssociatorAtVertexCMSTopTagCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('cmsTopTagCHS'), coneSize = 0.8)
process.secondaryVertexTagInfosCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCMSTopTagCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCMSTopTagCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCMSTopTagCHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCMSTopTagCHS.jetTagComputer = cms.string('combinedSecondaryVertexCMSTopTagCHS')
process.patJetsCMSTopTagCHS.addTagInfos = True
process.patJetsCMSTopTagCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfos')
    )

addJetCollection(
    process,
    labelName = 'CMSTopTagCHSSubjets',
    jetSource = cms.InputTag('cmsTopTagCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchCMSTopTagCHSSubjets.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagCHSSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagCHSSubjets.matched = 'ca8GenJets'#slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHSSubjets.matched = importantgenparticles

process.patJetsCMSTopTagCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCMSTopTagCHS" ),
    subjetSrc=cms.InputTag("patJetsCMSTopTagCHSSubjets")
      )


process.patMETs.addGenMET = False # There's no point in recalculating this, and we can't remake it since we don't have genParticles beyond |eta|=5

##########################################################################################
# PRODUCER

process.ttbsmAna = cms.EDFilter('TTBSMProducer',
                                wTagSrc = cms.InputTag('patJetsCA8CHSprunedPacked'),
                                topTagSrc = cms.InputTag('patJetsCMSTopTagCHSPacked'),
                                trigSrc = cms.InputTag('hltTriggerSummaryAOD'),
                                rhoSrc = cms.InputTag('kt6PFJets','rho'),#'fixedGridRhoAll','')#ak4PFJets', 'rho'),#kt6PFJets
                                genJetsSrc = cms.InputTag('ca8GenJets'),
				                        pvSrc = cms.InputTag(primaryvertices),
                                readTrig = cms.bool(False),
                                trigs = cms.vstring(
                                    ''
                                    ),
                                topTagParams = caTopTagParams.clone(
                                    tagName = cms.string('CATop')
                                    ),
                                wTagParams = boostedTopWTagParams.clone(
                                    yCut = cms.double(0.0)
                                    ),
                                jetScale = cms.double(0.0),   #
                                jetPtSmear = cms.double(0.1), # note these three are fractional
                                jetEtaSmear = cms.double(0.1),#
                                jecPayloads = cms.vstring([
                                    'START53_L1FastJet_AK7PFchs.txt',
                                    'START53_L2Relative_AK7PFchs.txt',
                                    'START53_L3Absolute_AK7PFchs.txt',
                                    #'START53_L2L3Residual_AK7PFchs.txt',
                                    'START53_Uncertainty_AK7PFchs.txt'
                                    ]),
                                pdfSet = cms.string("")

)


process.ttbsmAnaScaleDown = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetScale = cms.double(-0.05)
    )
process.ttbsmAnaScaleUp = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetScale = cms.double(0.05)
    )
process.ttbsmAnaPtSmearDown = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetPtSmear = cms.double(0.0)
    )
process.ttbsmAnaPtSmearUp = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetPtSmear = cms.double(0.2)
    )
process.ttbsmAnaEtaSmearDown = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetEtaSmear = cms.double(0.0)
    )
process.ttbsmAnaEtaSmearUp = process.ttbsmAna.clone(
    pdfSet = cms.string(""),
    jetEtaSmear = cms.double(0.2)
    )

##########################################################################################
# PDF

# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
      # Fix POWHEG if buggy (this PDF set will also appear on output, 
      # so only two more PDF sets can be added in PdfSetNames if not "")
      FixPOWHEG = cms.untracked.string(""),
      GenTag = cms.untracked.InputTag("prunedGenParticles"),
      PdfInfoTag = cms.untracked.InputTag("generator"),
      PdfSetNames = cms.untracked.vstring( "cteq66.LHgrid"
                                          ,"MRST2006nnlo.LHgrid"
                                           )
)

##########################################################################################
# PATHS


print 'Making the path'

process.p = cms.Path(
    #process.patTriggerDefaultSequence*
    #process.hltHighLevel
    #   process.ak4PFJets
    # *process.ca8PFJets
    # *process.ak15PFJets
    # process.ca8GenJets
    # *process.ca8GenJetsPruned
    # *process.chs
    # *process.ca8CHSJets
    #  *process.ca8CHSJetsPruned
    # *process.cmsTopTagCHS
    #*process.patJetsCA8CHS
    # *process.patJetCorrFactorsCA8CHS
    # process.patJetsCA8CHSpruned
    # *process.patJetsCA8CHSprunedSubjets
    # *process.patJetsCA8CHSprunedPacked
    # *process.patJetsCMSTopTagCHS
    # *process.patJetsCMSTopTagCHSSubjets
    # *process.patJetsCMSTopTagCHSPacked
    # *process.ttbsmAna
    # *process.ttbsmAnaScaleDown
    # *process.ttbsmAnaScaleUp
    # *process.ttbsmAnaPtSmearDown
    # *process.ttbsmAnaPtSmearUp
    # *process.ttbsmAnaEtaSmearDown
    # *process.ttbsmAnaEtaSmearUp
    #*process.pdfWeights
    )


process.outpath = cms.EndPath(process.out)

process.out.dropMetaData = cms.untracked.string("DROPPED")
