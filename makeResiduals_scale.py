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
ROOT.gSystem.Load("GeneralizedEndpoint_cc.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

def getAngle(v1,v2):
    vector1 = ROOT.TVector3(v1.Px(),v1.Py(),v1.Pz())
    vector2 = ROOT.TVector3(v2.Px(),v2.Py(),v2.Pz())
    
    return vector1.Angle(vector2)    

def whereIsDiMuon(eta1,eta2):
    abseta1 = abs(eta1)
    abseta2 = abs(eta2)
    
    if (abseta1 < 0.8):
        if (abseta2 < 0.8):    return "BB"
        elif (abseta2 < 1.2):  return "BO"
        elif (abseta2 < 2.45): return "BE"
    elif (abseta1 < 1.2):   
        if (abseta2 < 0.8):    return "BO"
        elif (abseta2 < 1.2):  return "OO"
        elif (abseta2 < 2.45): return "OE"
    elif (abseta1 < 2.45):  
        if (abseta2 < 0.8):    return "BE"
        elif (abseta2 < 1.2):  return "OE"
        elif (abseta2 < 2.45): return "EE"
    return "UU"

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

def correctMuonPt(oldmu):
    mu = oldmu
    mupt = ROOT.GeneralizedEndpoint().GeneralizedEndpointPt(mu.pt(),mu.charge(),mu.eta(),math.degrees(mu.phi()),0,1)

#    p4 = ROOT.Math.LorentzVector
    p4 = ROOT.reco.Candidate.PolarLorentzVector(mupt,mu.eta(),mu.phi(),0.10566)
    mu.setP4(p4)
    
    return mu


def makeResiduals(inputfilename, outputfilename):
    # create residual histograms
    DimuonMassGEN            = ROOT.TH1F("DimuonMassGEN"           ,"Dimu mass GEN", 600,  0,  6000)
    DimuonMassRECO           = ROOT.TH1F("DimuonMassRECO"          ,"Dimu mass RECO",600,  0,  6000)
    DimuonMassRes_0to500     = ROOT.TH1F("DimuonMassRes_0to500"    ,"Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_500to1000  = ROOT.TH1F("DimuonMassRes_500to1000" ,"Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_1000to1500 = ROOT.TH1F("DimuonMassRes_1000to1500","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_1500to2000 = ROOT.TH1F("DimuonMassRes_1500to2000","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_2000to2500 = ROOT.TH1F("DimuonMassRes_2000to2500","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_2500to3000 = ROOT.TH1F("DimuonMassRes_2500to3000","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_3000to3500 = ROOT.TH1F("DimuonMassRes_3000to3500","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_3500to4000 = ROOT.TH1F("DimuonMassRes_3500to4000","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_4000to4500 = ROOT.TH1F("DimuonMassRes_4000to4500","Dimu mass res", 100, -0.5, 0.5)
    DimuonMassRes_4500to5000 = ROOT.TH1F("DimuonMassRes_4500to5000","Dimu mass res", 100, -0.5, 0.5)
    
    DimuonMassRes_0to500By     = []
    DimuonMassRes_500to1000By  = []
    DimuonMassRes_1000to1500By = []
    DimuonMassRes_1500to2000By = []
    DimuonMassRes_2000to2500By = []
    DimuonMassRes_2500to3000By = []
    DimuonMassRes_3000to3500By = []
    DimuonMassRes_3500to4000By = []
    DimuonMassRes_4000to4500By = []
    DimuonMassRes_4500to5000By = []
    
    where_names = ["BB","BO","BE","OO","OE","EE","UU"]
    for n in where_names:
        DimuonMassRes_0to500By    .append( ROOT.TH1F(("DimuonMassRes_0to500"    +n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_500to1000By .append( ROOT.TH1F(("DimuonMassRes_500to1000" +n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_1000to1500By.append( ROOT.TH1F(("DimuonMassRes_1000to1500"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_1500to2000By.append( ROOT.TH1F(("DimuonMassRes_1500to2000"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_2000to2500By.append( ROOT.TH1F(("DimuonMassRes_2000to2500"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_2500to3000By.append( ROOT.TH1F(("DimuonMassRes_2500to3000"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_3000to3500By.append( ROOT.TH1F(("DimuonMassRes_3000to3500"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_3500to4000By.append( ROOT.TH1F(("DimuonMassRes_3500to4000"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_4000to4500By.append( ROOT.TH1F(("DimuonMassRes_4000to4500"+n),"Dimu mass res", 100, -0.5, 0.5))
        DimuonMassRes_4500to5000By.append( ROOT.TH1F(("DimuonMassRes_4500to5000"+n),"Dimu mass res", 100, -0.5, 0.5))
        
    
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

    for iev,event in enumerate(events):
#        if iev > 10: break
#        print "Event", iev 

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
##        for j,to in enumerate(triggerObjects.product()):
##            to.unpackPathNames(names);
##            print "Trigger object pt %6.2f eta %+5.3f phi %+5.3f  " % (to.pt(),to.eta(),to.phi())
##            print "         collection: ", to.collection()
##            print "         type ids: ", ", ".join([str(f) for f in to.filterIds()])
##            print "         filters: ", ", ".join([str(f) for f in to.filterLabels()])
##            pathslast = set(to.pathNames(True))
##            print "         paths:   ", ", ".join([("%s*" if f in pathslast else "%s")%f for f in to.pathNames()]) 
        
        # Gen Info
        selectedGenMuonsSt1 = []
        selectedGenMuonsSt23 = []
        selectedGenMuonsHP = []
        gen_particles = genparticles.product()
#        packed = packedgenparticles.product()

        for genp in gen_particles:
            if (abs(genp.pdgId())!=13): continue
            if (genp.status() == 23):
                selectedGenMuonsSt23.append(genp)
            if (genp.status() == 1):
                selectedGenMuonsSt1.append(genp)                
#                print 'saving muon with st1 pt: %f' %(genp.pt()) 
            if (genp.isHardProcess()):  ##ONLY STORE MUONS COMING FROM THE HARD PROCESS NO MATTER THE STATUS...
#                print 'saving muon with pt: %f' %(genp.pt()) 
                selectedGenMuonsHP.append(genp)                
                        
        selectedGenMuons = selectedGenMuonsHP
#        selectedGenMuons = selectedGenMuonsSt23
#        selectedGenMuons = selectedGenMuonsSt1

        if (len(selectedGenMuons) < 2): continue
        genMu1 = selectedGenMuons[0].p4()
        genMu2 = selectedGenMuons[1].p4()

        # Vertices
        if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4:
            #print "Event %8d has no good primary vertex." %(iev)
            continue
        else:
            PV = vertices.product()[0]
        
        # Muons
        selectedMuons = []
#        selectedMuonsB4 = []
        for j,muproxy in enumerate(muons.product()): 
            if abs(muproxy.eta()) > 2.4: continue
            mu = correctMuonPt(muproxy)
            
            if mu.pt() < 53: continue
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
        
        mu1 = selectedMuons[0].p4()
        mu2 = selectedMuons[1].p4()
        if (selectedMuons[0].charge()*selectedMuons[1].charge() > 0): continue
        if (getAngle(mu1,mu2)>3.1216): continue
        
        mass = (mu1+mu2).M()
  #      massB4 = (selectedMuonsB4[0].p4()+selectedMuonsB4[1].p4()).M()
        gen_mass = (genMu1+genMu2).M()
 #       print "Gen Mass = %f  ----> Mass before: %f, after: %f" %(gen_mass,massB4,mass)
        rdil = (mass / gen_mass) - 1

        DimuonMassGEN.Fill(gen_mass)
        DimuonMassRECO.Fill(mass)

        w = where_names.index(whereIsDiMuon(mu1.Eta(), mu2.Eta()))
        if   (gen_mass > 0    and gen_mass <= 500 ):  
            DimuonMassRes_0to500     .Fill(rdil)  
            DimuonMassRes_0to500By[w].Fill(rdil)
        elif (gen_mass > 500  and gen_mass <= 1000):  
            DimuonMassRes_500to1000 .Fill(rdil)   
            DimuonMassRes_500to1000By [w].Fill(rdil) 
        elif (gen_mass > 1000 and gen_mass <= 1500):  
            DimuonMassRes_1000to1500.Fill(rdil)   
            DimuonMassRes_1000to1500By[w].Fill(rdil) 
        elif (gen_mass > 1500 and gen_mass <= 2000):  
            DimuonMassRes_1500to2000.Fill(rdil)   
            DimuonMassRes_1500to2000By[w].Fill(rdil) 
        elif (gen_mass > 2000 and gen_mass <= 2500):  
            DimuonMassRes_2000to2500.Fill(rdil)   
            DimuonMassRes_2000to2500By[w].Fill(rdil) 
        elif (gen_mass > 2500 and gen_mass <= 3000):  
            DimuonMassRes_2500to3000.Fill(rdil)   
            DimuonMassRes_2500to3000By[w].Fill(rdil) 
        elif (gen_mass > 3000 and gen_mass <= 3500):  
            DimuonMassRes_3000to3500.Fill(rdil)   
            DimuonMassRes_3000to3500By[w].Fill(rdil) 
        elif (gen_mass > 3500 and gen_mass <= 4000):  
            DimuonMassRes_3500to4000.Fill(rdil)   
            DimuonMassRes_3500to4000By[w].Fill(rdil) 
        elif (gen_mass > 4000 and gen_mass <= 4500):  
            DimuonMassRes_4000to4500.Fill(rdil)   
            DimuonMassRes_4000to4500By[w].Fill(rdil) 
        elif (gen_mass > 4500 and gen_mass <= 5000):  
            DimuonMassRes_4500to5000.Fill(rdil)   
            DimuonMassRes_4500to5000By[w].Fill(rdil) 

        # printing out some info
        if iev > 0 and iev % 1000 == 0:
            print "Processed entry %8d of this tree" % (iev)
            
            
    # Now save the output to some random TFile and have a look later...
    print "Saving histograms"
        
    f = ROOT.TFile(outputfilename,"RECREATE");
    DimuonMassGEN .Write()
    DimuonMassRECO.Write()
    DimuonMassRes_0to500    .Write()
    DimuonMassRes_500to1000 .Write()
    DimuonMassRes_1000to1500.Write()
    DimuonMassRes_1500to2000.Write()
    DimuonMassRes_2000to2500.Write()
    DimuonMassRes_2500to3000.Write()
    DimuonMassRes_3000to3500.Write()
    DimuonMassRes_3500to4000.Write()
    DimuonMassRes_4000to4500.Write()
    DimuonMassRes_4500to5000.Write()
    for i,n in enumerate(where_names):
        DimuonMassRes_0to500By[i]    .Write()
        DimuonMassRes_500to1000By[i] .Write()
        DimuonMassRes_1000to1500By[i].Write()
        DimuonMassRes_1500to2000By[i].Write()
        DimuonMassRes_2000to2500By[i].Write()
        DimuonMassRes_2500to3000By[i].Write()
        DimuonMassRes_3000to3500By[i].Write()
        DimuonMassRes_3500to4000By[i].Write()
        DimuonMassRes_4000to4500By[i].Write()
        DimuonMassRes_4500to5000By[i].Write()
        
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
