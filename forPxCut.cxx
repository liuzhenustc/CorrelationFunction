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
TH1D *hEvent;
TH2D *hCen9VsPxL;
TH2D *hCen9VsPxR;
TH2D *hCen16VsPxL;
TH2D *hCen16VsPxR;

//----define-function---
void bookHistograms();
void writeHistograms(char* outFile);
Bool_t passPion(miniDst* evt);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2, Float_t *parL1, Float_t *parL2, Float_t *parR1, Float_t *parR2);

const Bool_t debug = 0;
//---global--var---------------
Float_t  vz = 0.;
Float_t  PxL = 0.;
Float_t  PxR = 0.;
Int_t    NTrks = 0;
Short_t  centrality9 = 0.;
Short_t  centrality16 = 0.;
Int_t    iran = 0;
Double_t reWeight = 1.;
//------------------------------

LorentzVec PiCandidate;
FloatVec vecPxL;
FloatVec vecPxR;

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
                cout<<"something wrong: "<<filename<<endl;
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

        hEvent->Fill(0.5);

        PiCandidate.clear();
        vecPxL.clear();
        vecPxR.clear();

        //--parameters--global-----
        NTrks = event->mNTrk;
        vz = event -> mTpcVz ;
        centrality9 = event -> mCentrality9;
        centrality16 = event -> mCentrality16;
        //----------------------------

        if(!passPion(event)) continue;//Get PiCandidate
        Int_t nPi = PiCandidate.size();
        if(nPi<2) continue;
        if(debug) cout<<"nPi: "<< nPi <<endl;
        hEvent->Fill(1.5);

        for(Int_t i=1; i<nPi; i++){
            comparePt(&PiCandidate[0],&PiCandidate[i],&vecPxL[0],&vecPxL[i],&vecPxR[0],&vecPxR[i]);
        }

        hCen9VsPxL->Fill(centrality9,vecPxL[0]);
        hCen9VsPxR->Fill(centrality9,vecPxR[0]);
        hCen16VsPxL->Fill(centrality16,vecPxL[0]);
        hCen16VsPxR->Fill(centrality16,vecPxR[0]);
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
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2, Float_t *parL1, Float_t *parL2, Float_t *parR1, Float_t *parR2)
{
    TLorentzVector pi1=*vec1;
    TLorentzVector pi2=*vec2;
    Float_t pt1 = pi1.Pt();
    Float_t pt2 = pi2.Pt();
    TLorentzVector vec;
    Float_t pxL;
    Float_t pxR;
    if(pt1<pt2) {
        vec=*vec1;
        *vec1=*vec2;
        *vec2=vec;
        pxL=*parL1;
        *parL1=*parL2;
        *parL2=pxL;
        pxR=*parR1;
        *parR1=*parR2;
        *parR2=pxR;
    }
}
//-------------------------------------------
Bool_t passPion(miniDst* evt)
{
    for(Int_t i=0;i<NTrks;i++){
        Float_t nSigmaPi = evt -> mnSigmaPi[i];
        if(!(nSigmaPi>0. && nSigmaPi<3.)) continue;
        Float_t pxL = evt -> mPxL[i];
        vecPxL.push_back(pxL);
        Float_t pxR = evt -> mPxR[i];
        vecPxR.push_back(pxR);
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
void bookHistograms()
{
    hEvent = new TH1D("hEvent","# of Event",5,0,5);
    hEvent->GetXaxis()->SetBinLabel(1,"all dimuon event");
    hEvent->GetXaxis()->SetBinLabel(2,"pass 2pion cut");
   
    hCen9VsPxL = new TH2D("hCen9VsPxL","PxL vs 9 Centrality Bin; Centrality; PxL",9,0,9,180,-70,20);
    hCen9VsPxR = new TH2D("hCen9VsPxR","PxR vs 9 Centrality Bin; Centrality; PxR",9,0,9,180,-70,20);
    hCen16VsPxL = new TH2D("hCen16VsPxL","PxL vs 16 Centrality Bin; Centrality; PxL",16,0,16,180,-70,20);
    hCen16VsPxR = new TH2D("hCen16VsPxR","PxR vs 16 Centrality Bin; Centrality; PxR",16,0,16,180,-70,20);
}
//--------------------------------------------
void writeHistograms(char* outFile)
{
    char buf[1024];
    sprintf(buf,"%s.histo.root",outFile);
    cout<<"Writing histograms into "<<buf<<endl;
    TFile *mFile = new TFile(buf,"recreate");
    mFile->cd();

    hEvent    -> Write(); 
    hCen9VsPxL -> Write(); 
    hCen9VsPxR -> Write(); 
    hCen16VsPxL -> Write(); 
    hCen16VsPxR -> Write(); 
}
