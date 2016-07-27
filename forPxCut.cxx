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

        hEvent->Fill(0.5);

        PiCandidate.clear();

        //--parameters--global-----
        NTrks = event->mNTrk;
        vz = event -> mTpcVz ;
        centrality9 = event -> mCentrality9;
        centrality16 = event -> mCentrality16;
        PxL = event -> mPxL;
        PxR = event -> mPxR;
        //----------------------------

        if(!passPion(event)) continue;//Get PiCandidate
        Int_t nPi = PiCandidate.size();
        if(nPi<2) continue;
        if(debug) cout<<"nPi: "<< nPi <<endl;
        hEvent->Fill(1.5);

        hCen9VsPxL->Fill(centrality9,PxL);
        hCen9VsPxR->Fill(centrality9,PxR);
        hCen16VsPxL->Fill(centrality16,PxL);
        hCen16VsPxR->Fill(centrality16,PxR);
    }

    timer.Stop();

    writeHistograms(outFile);
    delete chain;
    cout<<"end of program"<<endl;
    timer.Print();
    return 0;
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
