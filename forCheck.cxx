#include <cstdio> 
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "sys/types.h"
#include "dirent.h"
#include "math.h"
#include "string.h"

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "miniDst.h"
#include "THnSparse.h"
using namespace std;
#endif

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
//typedef vector<Bool_t> BoolVec;
//typedef vector<Int_t> IntVec;
typedef vector<Float_t> FloatVec;
//typedef vector<Double_t> DoubleVec;
typedef vector<TLorentzVector> LorentzVec;
#else
//typedef vector<Bool_t, allocator<Bool_t>> BoolVec;
//typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<Float_t, allocator<Float_t>> FloatVec;
//typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
#endif

const Float_t Mpi = 0.13957;//GeV

//--histogram
TH2D *hVertexYvsX;
TH2D *hVpdVzvsTpcVz;
TH1D *hVzDiff;
TH1D *hRefMult;
TH2D *hGRefMultvsGRefMultCorr;
TH1D *hCentrality;

TH1D *hnPiTrig;
TH1D *hTrigPt;
TH1D *hAssocPt;

THnSparseD *hDetaDphiVsPtReal;
THnSparseD *hDetaDphiVsPtMix;

const Int_t dimCor = 4;
const Int_t nCorBins[dimCor]={400,400,32,96};//trigPt;AssocPt;dEta;dPhi
const Double_t CorlowBins[dimCor]={0,0,-1.6,-0.5*3.1415926};
const Double_t CorupBins[dimCor]={20,20,1.6,1.5*3.1415926};

//----define-function---
void bookHistograms();
void writeHistograms(char* outFile);
Bool_t passEvent(miniDst* evt);
Bool_t passPion(miniDst* evt);
Bool_t passPionForReal(LorentzVec vec);
Bool_t IsPxEvent(Float_t Px,Short_t cen,Float_t PxCut);
Bool_t GetPtPointer(LorentzVec vecCandidate);
Bool_t makeMixPairs(LorentzVec vecTrig,THnSparseD *hn);
void copyToBuffer(LorentzVec vecCandidate);
Bool_t makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void copyToBuffer(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2);

const Bool_t debug = 0;
//---global--var---------------
Float_t  vz = 0.;//tpc vz
Int_t    NTrks = 0;
Short_t  centrality = 0.;
Int_t    iran = 0;
Double_t reWeight = 1.;
//------------------------------

LorentzVec PiCandidate;
LorentzVec PiTrig;
LorentzVec PiAssoc;

//------var-used-for-mixevent-------
const Int_t mCenBins = 9; //default:16
const Int_t mVzBins = 20; //default:20
const Int_t mPtBins = 8; //1-2;2-3;3-4;4-5;5-6,6-7,7-8,>8
const Int_t mMaxEvtsInBuffer = 100; //default:50
const Double_t VzCut = 100.;

Int_t cenPointer, vzPointer,ptPointer;
Int_t nEvtsInBuffer[mCenBins][mVzBins][mPtBins];
Bool_t bufferFullFlag[mCenBins][mVzBins][mPtBins];
LorentzVec bufferPi[mCenBins][mVzBins][mPtBins][mMaxEvtsInBuffer];
//----------------------------------

int main(int argc, char** argv)
{
    if(argc!=1&&argc!=3) return -1;

    TString inFile="test.list";
    char outFile[1024];
    sprintf(outFile,"test/test");
    if(argc==3){
        inFile = argv[1];
        sprintf(outFile,"%s",argv[2]);
    }

    //---open files--and--add to chain---//
    TChain *chain = new TChain("miniDst");

    Int_t ifile=0;
    char filename[512];
    ifstream *inputStream = new ifstream;
    inputStream->open(inFile.Data());
    if (!(inputStream)) {
        printf("can not open list file\n");
        return 0;
    }
    for(;inputStream->good();){
        inputStream->getline(filename,512);
        if(inputStream->good()) {
            TFile *ftmp = new TFile(filename);
            if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
                cout<<"something wrong"<<endl;
            } else {
                if(debug) cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
                chain->Add(filename);
                ifile++;
            }
            delete ftmp;
        }
    }
    delete inputStream;

    bookHistograms();

    miniDst *event = new miniDst(chain);

    Int_t nEvts = chain->GetEntries();
    cout<<nEvts<<" events"<<endl;
    int perevt = nEvts/100;

    //-----timer---------
    TStopwatch timer;
    timer.Start();

    //----loop--events---//
    for(Int_t i=0;i<nEvts;i++) {
        if(i%perevt ==0) cout<<"Looking at "<<((float)i/nEvts)*100<<" % .."<<endl;
        event->GetEntry(i);

        PiTrig.clear();
        PiAssoc.clear();
        PiCandidate.clear();

        //--parameters--global-----
        NTrks = event->mNTrk;
        vz = event -> mTpcVz ;
        centrality = event -> mCentrality9;
        //----------------------------

        if(!passEvent(event)) continue;
        cenPointer = cenPointer;
        vzPointer = (Int_t)( (vz+VzCut)/(2*VzCut)*mVzBins );
        if(debug) cout<<"centrality: "<<centrality<<"  cenpointer: "<<cenPointer<<"  vzPointer:  "<<vzPointer<<endl;

        if(!passPion(event)) continue;//Get PiCandidate
        Int_t nPi = PiCandidate.size();
        if(nPi<2) continue;
        if(debug) cout<<"nPi: "<< nPi <<endl;

        //GetPtPointer(PiCandidate);
        //if(debug) cout<<"ptPointer: "<<ptPointer<<endl;

        if(!passPionForReal(PiCandidate)) continue;//Get Trigger&Assoc

        hTrigPt->Fill(PiCandidate[0].Pt());
        for(unsigned int i=1;i<PiCandidate.size();i++) hAssocPt->Fill(PiCandidate[i].Pt());//associate
        hnPiTrig->Fill(7.5);//trigger

        if(debug) cout<<"# of Trigger Particle "<<PiTrig.size()<<endl;
        if(debug) cout<<"# of Assoc Particle "<<PiAssoc.size()<<endl;

        if(!makeRealPairs(PiTrig,PiAssoc,hDetaDphiVsPtReal)) continue; 
    }

    timer.Stop();

    writeHistograms(outFile);
    delete chain;
    cout<<"end of program"<<endl;
    timer.Print();
    return 0;
}
//------------------------------------------
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2)
{
    TLorentzVector pi1=*vec1;
    TLorentzVector pi2=*vec2;
    Float_t pt1 = pi1.Pt();
    Float_t pt2 = pi2.Pt();
    TLorentzVector vec;
    if(pt1<pt2) {
        vec=*vec1;
        *vec1=*vec2;
        *vec2=vec;
    }

}
//------------------------------------------
Bool_t passEvent(miniDst* evt)
{
    Double_t vpdVz = evt->mVpdVz;
    Double_t vzDiff = vz - vpdVz;
    Double_t grefMult = evt->mGRefMult;
    Double_t grefMultCorr = evt->mGRefMultCorr;

    hVpdVzvsTpcVz->Fill(vz,vpdVz);
    hVzDiff->Fill(vzDiff);
    hGRefMultvsGRefMultCorr->Fill(grefMultCorr,grefMult);
    hCentrality->Fill(centrality);

    return kTRUE;
}
//-------------------------------------------
Bool_t passPion(miniDst* evt)
{
    for(Int_t i=0;i<NTrks;i++){
        Float_t nSigmaPi = evt -> mnSigmaPi[i];
        if(!(nSigmaPi>0. && nSigmaPi<3.)) continue;
        Float_t pt = evt -> mTrkPt[i];
        Float_t eta = evt -> mTrkEta[i];
        Float_t phi = evt -> mTrkPhi[i];
        TLorentzVector fourmom(0,0,0,0);
        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
        PiCandidate.push_back(fourmom);
    }
    if(debug)cout<<"passPion PiCandidate: "<<PiCandidate.size()<<endl;
    return kTRUE;
}
//-------------------------------------------
Bool_t passPionForReal(LorentzVec vec)
{
    Int_t nvec = vec.size();
    if(nvec<2) {
        if(debug) cout<<"PiCandidate error nPiCandidate: "<<nvec<<endl;
        return kFALSE;
    }
    for(Int_t i=1; i<nvec; i++){
        comparePt(&vec[0],&vec[i]);
    }

    PiTrig.push_back(vec[0]);
    if(debug)cout<<"passPionForReal Pass Trigger: "<<PiTrig[0].Pt()<<endl;
    
    for(Int_t i=1;i<nvec ; i++){
        PiAssoc.push_back(vec[i]);
    }
    return kTRUE;
}
//-------------------------------------------
Bool_t IsPxEvent(Float_t Px,Short_t cen,Float_t PxCut)
{
    if(Px < PxCut) return kTRUE;
    else return kFALSE;
}
//-------------------------------------------
Bool_t makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn)
{
    Int_t nAssoc = vecAssoc.size();
    if(nAssoc==0) return kFALSE;

    Float_t TrigPt = vecTrig[0].Pt();
    Float_t TrigEta = vecTrig[0].Eta();
    Float_t TrigPhi = vecTrig[0].Phi();

    for(Int_t i =0;i<nAssoc;i++){
        Float_t AssocPt = vecAssoc[i].Pt();
        Float_t AssocEta = vecAssoc[i].Eta();
        Float_t AssocPhi = vecAssoc[i].Phi();

        Float_t dEta = AssocEta - TrigEta;
        Float_t dPhi = AssocPhi - TrigPhi;
        if(dPhi < -TMath::Pi()/2.) dPhi += 2*TMath::Pi();
        if(dPhi > 3*TMath::Pi()/2.) dPhi -= 2*TMath::Pi();

        Double_t fill[]={TrigPt,AssocPt,dEta,dPhi};
        hn->Fill(fill);
    }
    return kTRUE;
}
//-------------------------------------------
Bool_t GetPtPointer(LorentzVec vecCandidate)
{
    Int_t nVec = vecCandidate.size();
    for(Int_t i=0;i<nVec;i++){
        comparePt(&vecCandidate[0],&vecCandidate[i]);
    }
    Float_t candidatePt = vecCandidate[0].Pt();
    if(1.<candidatePt && candidatePt<2.) ptPointer=0;
    else if(2.<candidatePt && candidatePt<3.) ptPointer=1;
    else if(3.<candidatePt && candidatePt<4.) ptPointer=2;
    else if(4.<candidatePt && candidatePt<5.) ptPointer=3;
    else if(5.<candidatePt && candidatePt<6.) ptPointer=4;
    else if(6.<candidatePt && candidatePt<7.) ptPointer=5;
    else if(7.<candidatePt && candidatePt<8.) ptPointer=6;
    else if(8.<candidatePt) ptPointer=7;

    if(debug)cout<<"GetPtPointer: "<<candidatePt<<"  prPointer: "<<ptPointer<<endl;

    return kFALSE;
} 
//-------------------------------------------
Bool_t makeMixPairs(LorentzVec vecTrig,THnSparseD *hn)
{
    Int_t nRealtrig = vecTrig.size();
    if(nRealtrig!=1){
        cout<<"Trigger Particle error"<<endl;
        return kFALSE;
    }

    Float_t TrigPt = vecTrig[0].Pt();
    Float_t TrigEta = vecTrig[0].Eta();
    Float_t TrigPhi = vecTrig[0].Phi();

    for(Int_t k=0;k<mPtBins;k++){
        ptPointer = k;
        for(Int_t evtPointer =0;evtPointer < nEvtsInBuffer[cenPointer][vzPointer][ptPointer];evtPointer++){
            Int_t nbufferPi = bufferPi[cenPointer][vzPointer][ptPointer][evtPointer].size();
            for(Int_t i=0;i<nbufferPi;i++){
                Float_t BufferPt = bufferPi[cenPointer][vzPointer][ptPointer][evtPointer][i].Pt();
                Float_t BufferEta = bufferPi[cenPointer][vzPointer][ptPointer][evtPointer][i].Eta();
                Float_t BufferPhi = bufferPi[cenPointer][vzPointer][ptPointer][evtPointer][i].Phi();

                Float_t dEta = BufferEta - TrigEta;
                Float_t dPhi = BufferPhi - TrigPhi;
                if(dPhi < -TMath::Pi()/2.) dPhi += 2*TMath::Pi();
                if(dPhi > 3*TMath::Pi()/2.) dPhi -= 2*TMath::Pi();
                Double_t fill[]={TrigPt,BufferPt,dEta,dPhi};
                hDetaDphiVsPtMix->Fill(fill);
            }
        }
    }
    return kTRUE;
}
//-------------------------------------------
void copyToBuffer(LorentzVec vecCandidate)
{
    if(nEvtsInBuffer[cenPointer][vzPointer][ptPointer]>=mMaxEvtsInBuffer){
        bufferFullFlag[cenPointer][vzPointer][ptPointer] = kTRUE;
    }//full flag
    TRandom3 *gRandom = new TRandom3(iran++);
    Int_t evtPointer = -1;
    if(bufferFullFlag[cenPointer][vzPointer][ptPointer]) evtPointer = (Int_t)gRandom->Uniform(0,mMaxEvtsInBuffer-1.e-6);
    else evtPointer = nEvtsInBuffer[cenPointer][vzPointer][ptPointer];
    delete gRandom;

    bufferPi[cenPointer][vzPointer][ptPointer][evtPointer].clear();
    Int_t nVec = vecCandidate.size();
    if(debug) cout<<" evtPointer: "<<evtPointer<<endl;
    for(Int_t i=0;i<nVec;i++){
        bufferPi[cenPointer][vzPointer][ptPointer][evtPointer].push_back(vecCandidate[i]);
    }
    if(nEvtsInBuffer[cenPointer][vzPointer][ptPointer]<mMaxEvtsInBuffer){
        nEvtsInBuffer[cenPointer][vzPointer][ptPointer]++;
    }

}
//-------------------------------------------
void bookHistograms()
{
    hVertexYvsX = new TH2D("hVertexYvsX","hVertexYvsX; VertexX (cm); VertexY (cm)",500,-2.5,2.5,500,-2.5,2.5);
    hVpdVzvsTpcVz = new TH2D("hVpdVzvsTpcVz","hVpdVzvsTpcVz; TPC VertexZ (cm); VPD VertexZ (cm)",2000,-200,200,2000,-200,200);
    hVzDiff = new TH1D("hVzDiff","hVzDiff;Vz_{TPC} - Vz_{VPD} (cm)",200,-10,10);
    hRefMult = new TH1D("hRefMult","hRefMult;refMult",1000,0,1000);
    hGRefMultvsGRefMultCorr = new TH2D("hGRefMultvsGRefMultCorr","hGRefMultvsGRefMultCorr;grefMultCorr;grefMult",1000,0,1000,1000,0,1000);
    hCentrality = new TH1D("hCentrality","hCentrality; mCentrality",16,0,16);

    hnPiTrig = new TH1D("hnPiTrig","# of Triggered pion",10,0,10);
    hnPiTrig->GetXaxis()->SetBinLabel(1,"nPxLeft");
    hnPiTrig->GetXaxis()->SetBinLabel(2,"nPxLeftClose");
    hnPiTrig->GetXaxis()->SetBinLabel(3,"nPxLeftFar");
    hnPiTrig->GetXaxis()->SetBinLabel(4,"nPxRight");
    hnPiTrig->GetXaxis()->SetBinLabel(5,"nPxRightClose");
    hnPiTrig->GetXaxis()->SetBinLabel(6,"nPxRightFar");
    hnPiTrig->GetXaxis()->SetBinLabel(7,"nBothLeftRight");
    hnPiTrig->GetXaxis()->SetBinLabel(8,"nEvents");

    hTrigPt = new TH1D("hTrigPt"," p_{T}^{trig} event level ;p_{T}^{trig} (GeV/c)",400,0,20);
    hAssocPt = new TH1D("hAssocPt","hAssocPt;p_{T}^{assoc} (GeV/c)",400,0,20);

    hDetaDphiVsPtReal = new THnSparseD("hDetaDphiVsPtReal","Realevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtMix = new THnSparseD("hDetaDphiVsPtMix","mixevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
}
//--------------------------------------------
void writeHistograms(char* outFile)
{
    char buf[1024];
    sprintf(buf,"%s.histo.root",outFile);
    cout<<"Writing histograms into "<<buf<<endl;
    TFile *mFile = new TFile(buf,"recreate");
    mFile->cd();

    hVertexYvsX->Write();
    hVpdVzvsTpcVz->Write();
    hVzDiff->Write();
    hRefMult->Write();
    hGRefMultvsGRefMultCorr->Write();
    hCentrality->Write();

    hnPiTrig->Write();
    hTrigPt->Write();
    hAssocPt->Write();

    hDetaDphiVsPtReal -> Write();
    hDetaDphiVsPtMix -> Write();
}
