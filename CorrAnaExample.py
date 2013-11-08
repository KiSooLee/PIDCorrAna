import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
                )
                            )
#import sys
#import os
#inputfilename = open("MBCorrelation.txt", "r")
#ffrom=int(sys.argv[2])
#fto=int(sys.argv[3])
#for j in range(1,ffrom+1):
#    inputfilename.readline()
#for i in range(ffrom,fto):
#    process.source.fileNames.append(inputfilename.readline())

process.demo = cms.EDAnalyzer('CorrAnaExample',
multMax = cms.untracked.double(30),
                multMin = cms.untracked.double(20)
                              )

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('yourfilename.root')
                                   )


process.p = cms.Path(process.demo)
