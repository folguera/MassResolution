#!/usr/bin/python

# import ROOT in batch mode
import sys,getopt

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
import math
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
#ROOT.gSystem.Load("GeneralizedEndpoint_cc.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

import json

with open('Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_MuonPhys.txt') as data_file:    
    jsonfile = json.load(data_file)


ptbins = [20,30,50,70,100,150,200,250,300,600]
categories = ["BB","BEm","BEp"]

def jsonContainsEvent (event,lumi):

   # if the jsonVec is empty, then no JSON file was provided so all
   # events should pass
   if not jsonfile: return True

   if str(event) in jsonfile: 
       for ls in jsonfile[str(event)]:
           if lumi >= ls[0] and lumi <= ls[1]: return True
   else: 
       return False
        
def getAngle(v1,v2):
    vector1 = ROOT.TVector3(v1.Px(),v1.Py(),v1.Pz())
    vector2 = ROOT.TVector3(v2.Px(),v2.Py(),v2.Pz())
    
    return vector1.Angle(vector2)    

def whereIsDiMuon(eta1,eta2):   
    if   (abs(eta1) <  1.2 and abs(eta2) <  1.2): return "BB"
    elif (    eta1  < -1.2 or      eta2  < -1.2): return "BEm"
    else:                                         return "BEp"

def getFiles(inputfile):
    if (".root" in inputfile): return inputfile
    else:
        files = []
        with open(inputfile, 'r') as f:
            for line in f:
                filename = "root://eoscms//eos/cms/"+line.rstrip('\n')
                files.append(filename)
        f.closed
        return files

def isAncestor(a,p) :
    if a == p : 
        return True
    for i in xrange(0,p.numberOfMothers()) :
        if isAncestor(a,p.mother(i)) :
            return True
    return False

def getP4withTrkPt(mu):
    if abs(mu.eta()) > 0.9: 
        mupt = mu.innerTrack().pt()
    else:
        mupt = mu.pt()
    
    p4 = ROOT.TLorentzVector(0,0,0,0)
    p4.SetPtEtaPhiM(mupt, mu.eta(),mu.phi(),0.10566)
####    p4 = ROOT.reco.Candidate.PolarLorentzVector(
    return p4

def makeResiduals(inputfilename, outputfilename):
    # create residual histograms
    DimuonMassRECO           = ROOT.TH1F("DimuonMassRECO"          ,"Dimu mass RECO",200,  0,  200)
    

    DiMuonMass = [[0 for x in range(len(ptbins)-1)] for y in range(len(categories))]
    PtErr  = [[0 for x in range(len(ptbins)-1)] for y in range(len(categories))]

    for i,pt in enumerate(ptbins):
        for c,cat in enumerate(categories):
            if pt==600: break
            DiMuonMass[c][i] = ROOT.TH1F("dimuonMass_%s_pt%sto%s" %(cat, ptbins[i],ptbins[i+1]),"dimuon mass res",80,60.,120.)
            PtErr     [c][i] = ROOT.TH1F("PtErr_%s_pt%sto%s" %(cat, ptbins[i],ptbins[i+1]),"pt error",100, 0., 0.3)

    # read information 
    muons, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
    vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlineSlimmedPrimaryVertices"
    verticesScore = Handle("edm::ValueMap<float>")                                                                                        
    genparticles, genLabel  = Handle ("std::vector<reco::GenParticle>"), "prunedGenParticles"
#    packedgenparticles, labelPacked = Handle ("std::vector<pat::PackedGenParticle>"), "packedGenParticles"

    triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
    triggerObjects, triggerObjectLabel  = Handle("std::vector<pat::TriggerObjectStandAlone>"), "selectedPatTrigger"
    triggerPrescales, triggerPrescaleLabel  = Handle("pat::PackedTriggerPrescales"), "patTrigger"

    # open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
    files = getFiles(inputfilename)
    events = Events( files )

    isData = False
    if('Run2016' in inputfilename): isData = True
    lumi = 0 
    run = 0
    checkLumi = False 
    for iev,event in enumerate(events):
#        if iev > 4000: break
#        print "Event", iev 
        if isData: 
            if lumi != event.eventAuxiliary().luminosityBlock(): 
                lumi = event.eventAuxiliary().luminosityBlock()
                checkLumi = True
            if run  != event.eventAuxiliary().run(): 
                run = event.eventAuxiliary().run()
                checkLumi = True

            if checkLumi:
                if not jsonContainsEvent(run, lumi):
                    continue

        event.getByLabel(muonLabel, muons)
        event.getByLabel(vertexLabel, vertices)
        event.getByLabel(vertexLabel, verticesScore)
        event.getByLabel(genLabel, genparticles)
        event.getByLabel(triggerBitLabel, triggerBits)
#        event.getByLabel(labelPacked, packedgenparticles)
#        event.getByLabel(triggerObjectLabel, triggerObjects)
#        event.getByLabel(triggerPrescaleLabel, triggerPrescales)
        
        names = event.object().triggerNames(triggerBits.product())
        passTrigger=False
        for i in xrange(triggerBits.product().size()):
            if ("HLT_Mu27_v" in names.triggerName(i)): passTrigger = passTrigger or triggerBits.product().accept(i) 
            if ("HLT_Mu50_v" in names.triggerName(i)): passTrigger = passTrigger or triggerBits.product().accept(i) 
            if ("HLT_Mu45_eta2p1_v" in names.triggerName(i)): passTrigger = passTrigger or triggerBits.product().accept(i) 

        if not (passTrigger): continue

        # Vertices
        if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
            #print "Event %8d has no good primary vertex." %(iev)
            continue
        else:
            PV = vertices.product()[0]
        
        # Muons
        selectedMuons = []
#        selectedMuonsB4 = []
        for j,mu in enumerate(muons.product()):
            if abs(mu.eta()) > 2.4: continue 
            if not mu.isGlobalMuon(): continue
            if not mu.isTrackerMuon(): continue
            if mu.globalTrack().hitPattern().numberOfValidMuonHits() == 0: continue
            if mu.globalTrack().hitPattern().numberOfValidPixelHits() == 0: continue
            if mu.globalTrack().hitPattern().trackerLayersWithMeasurement() <= 5: continue
            if mu.muonBestTrack().ptError()/mu.muonBestTrack().pt() > 0.3: continue
            if abs(mu.muonBestTrack().dxy(vertices.product()[0].position())) > 0.2: continue
            if mu.isolationR03().sumPt/mu.innerTrack().pt() > 0.1: continue
            
            selectedMuons.append(mu)
                    
        #apply dimuon selection
        if (len(selectedMuons) < 2): continue
        
        mu1 = getP4withTrkPt(selectedMuons[0])
        mu2 = getP4withTrkPt(selectedMuons[1])
        pt1 = mu1.Pt()
        pt2 = mu2.Pt()
        pterr1 = selectedMuons[0].muonBestTrack().ptError()/selectedMuons[0].muonBestTrack().pt()
        pterr2 = selectedMuons[1].muonBestTrack().ptError()/selectedMuons[1].muonBestTrack().pt()

        if (selectedMuons[0].charge()*selectedMuons[1].charge() > 0): continue
        if (getAngle(mu1,mu2)>3.1216): continue
        
        mass = (mu1+mu2).M()
        
        DimuonMassRECO.Fill(mass)
        
        cat = categories.index(whereIsDiMuon(mu1.Eta(),mu2.Eta()))
        for i,pt in enumerate(ptbins):
            if (pt==600):
                if (pt1 > ptbins[i] or pt2 > ptbins[i]): DiMuonMass[cat][i-1].Fill(mass)
#                    if cat=="BB" : DimuonMass_BB [i-1].Fill(mass)
#                    if cat=="BEp": DimuonMass_BEp[i-1].Fill(mass)
#                    if cat=="BEm": DimuonMass_BEm[i-1].Fill(mass)                                        
                if (pt1 > ptbins[i]):                    PtErr[cat][i-1].Fill(pterr1)
#                    if cat=="BB" : PtErr_BB [i-1].Fill(pterr1)
#                    if cat=="BEp": PtErr_BEp[i-1].Fill(pterr1)
#                    if cat=="BEm": PtErr_BEm[i-1].Fill(pterr1)                                        
                if (pt2 > ptbins[i]):                    PtErr[cat][i-1].Fill(pterr2)
#                    if cat=="BB" : PtErr_BB [i-1].Fill(pterr2)
#                    if cat=="BEp": PtErr_BEp[i-1].Fill(pterr2)
#                    if cat=="BEm": PtErr_BEm[i-1].Fill(pterr2)                                        
            else:
                if (pt1 > ptbins[i] and pt1 < ptbins[i+1]) or (pt2 > ptbins[i] and pt2 < ptbins[i+1]):  DiMuonMass[cat][i].Fill(mass)
                ## FOR PT
                if (pt1 > ptbins[i] and pt1 < ptbins[i+1]):  PtErr[cat][i].Fill(pterr1)
                if (pt2 > ptbins[i] and pt2 < ptbins[i+1]):  PtErr[cat][i].Fill(pterr2)
##                if (pt1 > ptbins[i] and pt1 < ptbins[i+1]) or (pt2 > ptbins[i] and pt2 < ptbins[i+1]):  
##                    if cat=="BB" : DimuonMass_BB [i].Fill(mass)
##                    if cat=="BEp": DimuonMass_BEp[i].Fill(mass)
##                    if cat=="BEm": DimuonMass_BEm[i].Fill(mass)                                        
##                if (pt1 > ptbins[i] and pt1 < ptbins[i+1]):  
##                    if cat=="BB" : PtErr_BB [i].Fill(pterr1)
##                    if cat=="BEp": PtErr_BEp[i].Fill(pterr1)
##                    if cat=="BEm": PtErr_BEm[i].Fill(pterr1)                                        
##                if (pt2 > ptbins[i] and pt2 < ptbins[i+1]):  
##                    if cat=="BB" : PtErr_BB [i].Fill(pterr2)
##                    if cat=="BEp": PtErr_BEp[i].Fill(pterr2)
##                    if cat=="BEm": PtErr_BEm[i].Fill(pterr2)                                        
                             
        # printing out some info
        if iev > 0 and iev % 1000 == 0:
            print "Processed entry %8d of this tree" % (iev)
            
            
    # Now save the output to some random TFile and have a look later...
    print "Saving histograms"        
    f = ROOT.TFile(outputfilename,"RECREATE");    
    DimuonMassRECO.Write()
    for i,pt in enumerate(ptbins):
        for c,cat in enumerate(categories):
            if pt==600: break
            DiMuonMass[c][i].Write()
            PtErr     [c][i].Write()                
    f.Close()

#### ========= MAIN =======================
def main(argv):
    inputfile = 'root://eoscms//eos/cms//store/mc/RunIIFall15MiniAODv2/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/247BB932-F7B8-E511-B2F3-D8D385AF891A.root'
    outputfile = 'histos.root'
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print 'makeResiduals.py -i <inputfile> -o <outputfile>'
            sys.exit(2)
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print "Running on: %s" %(inputfile)
    print "Saving result in: %s" %(outputfile)
    makeResiduals(inputfile,outputfile)        
        
if __name__ == "__main__":
    main(sys.argv[1:])
    print "DONE"
