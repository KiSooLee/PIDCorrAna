// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

#define PI 3.1416
using namespace std;

class CorrAnaExample : public edm::EDAnalyzer {
public:
    explicit CorrAnaExample(const edm::ParameterSet&);
    ~CorrAnaExample();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    TH1D* hMult_selected;
    TH1D* hMult;
    TH2D* hSignal;
    TH2D* hBackground;
    //TH2D* hCorrelation;
    
    vector<TVector3> *pVect_trg;
    vector< vector<TVector3> > *pVectVect_trg;
    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_trg_;
    double ptMax_trg_;
    double ptMin_ass_;
    double ptMax_ass_;
    int bkgFactor_;
    double multMax_;
    double multMin_;
    edm::Service<TFileService> fs;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CorrAnaExample::CorrAnaExample(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_trg_ = iConfig.getUntrackedParameter<double>("ptMin_trg", 0.3);
    ptMax_trg_ = iConfig.getUntrackedParameter<double>("ptMax_trg", 3.0);
    ptMin_ass_ = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_ = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 10);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
}


CorrAnaExample::~CorrAnaExample()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CorrAnaExample::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel("offlinePrimaryVertices",vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;
    
    pVect_trg = new vector<TVector3>;
    pVect_ass = new vector<TVector3>;
    //----- loop over tracks -----
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel("generalTracks",tracks);
    
    //track selection
    int nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult_selected->Fill(nMult_ass_good);
   
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double phi = trk.phi();
        double pt  = trk.pt();
        
        TVector3 pvector;
        pvector.SetPtEtaPhi(pt,eta,phi);
        if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptMax_trg_ && trk.pt()>=ptMin_trg_) pVect_trg->push_back(pvector);
        if(trk.eta()<=etaMax_ass_ && trk.eta()>=etaMin_ass_ && trk.pt()<=ptMax_ass_ && trk.pt()>=ptMin_ass_) pVect_ass->push_back(pvector);
    }
    // Calculating signal
    int nMult_trg = (int)pVect_trg->size();
    int nMult_ass = (int)pVect_ass->size();
    hMult->Fill(nMult_trg);
    
    for(int ntrg=0;ntrg<nMult_trg;ntrg++)
    {
        TVector3 pvector_trg = (*pVect_trg)[ntrg];
        double eta_trg = pvector_trg.Eta();
        double phi_trg = pvector_trg.Phi();
        
        for(int nass=0;nass<nMult_ass;nass++)
        {
            TVector3 pvector_ass = (*pVect_ass)[nass];
            double eta_ass = pvector_ass.Eta();
            double phi_ass = pvector_ass.Phi();
            
            double deltaEta=eta_ass-eta_trg;
            double deltaPhi=phi_ass-phi_trg;
            if(deltaPhi>PI)
                deltaPhi=deltaPhi-2*PI;
            if(deltaPhi<-PI)
                deltaPhi=deltaPhi+2*PI;
            if(deltaPhi>-PI && deltaPhi<-PI/2.)
                deltaPhi=deltaPhi+2*PI;
            
            //if(deltaEta==0 && deltaPhi==0) continue;
            if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue;
            hSignal->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
        }
    }
    pVectVect_trg->push_back(*pVect_trg);
    pVectVect_ass->push_back(*pVect_ass);
    zvtxVect->push_back(bestvz);
    
    delete pVect_trg;
    delete pVect_ass;
}
}

// ------------ method called once each job just before starting event loop  ------------
void
CorrAnaExample::beginJob(){
    
    TH1D::SetDefaultSumw2();
    
    hMult_selected = fs->make<TH1D>("mult_selected",";N",300,0,300);
    hMult = fs->make<TH1D>("mult",";N",1500,0,1500);
    hSignal = fs->make<TH2D>("signal",";#Delta#eta;#Delta#phi",
                             33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    hBackground = fs->make<TH2D>("background",";#Delta#eta;#Delta#phi",
                                 33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    //hCorrelation = fs->make<TH2D>("correlation",";#Delta#eta;#Delta#phi",
    //                              27,-4.05,4.05,33,-(0.5+1/32)*PI,(1.5+1/32)*PI);
    pVectVect_trg = new vector< vector<TVector3> >;
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
}

// ------------ method called once each job just after ending the event loop  ------------
void
CorrAnaExample::endJob() {
    // Calculate background
    int nevttotal_trg = (int)pVectVect_trg->size();
    int nevttotal_ass = (int)pVectVect_ass->size();
    
    for(int nround=0;nround<bkgFactor_;nround++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(nevttotal_trg);
            if(nevt_trg == nevt_ass) { nevt_ass--; continue; }
            if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                nevt_ass--;
                ncount++;
                if(ncount>5000) {nevt_ass++; ncount = 0;}
                continue; }
            
            vector<TVector3> pVectTmp_trg = (*pVectVect_trg)[nevt_trg];
            vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
            int nMult_trg = pVectTmp_trg.size();
            int nMult_ass = pVectTmp_ass.size();
            
            for(int ntrg=0;ntrg<nMult_trg;ntrg++)
            {
                TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                double phi_trg = pvectorTmp_trg.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta();
                    double phi_ass = pvectorTmp_ass.Phi();
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue;
                    
                    hBackground->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                }
            }
        }
    }
    
    /*int nEvent = hMult->Integral(3,10000);
    double Bz = hBackground->GetBinContent(275);
    hCorrelation->Add(hSignal);
    hCorrelation->Scale(Bz);
    hCorrelation->Divide(hBackground);
    hCorrelation->Scale(1.0/nEvent);*/
}

//define this as a plug-in
DEFINE_FWK_MODULE(CorrAnaExample);
