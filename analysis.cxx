//-----------------------------------------------------//
//distinguish PxL and PxR close and far region 
//include mix event
//-----------------------------------------------------//

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

TH1D *hEvent;
TH1D *hnPiTrig;

TH3D *hTrigPtEtaPhiPxLClose;
TH3D *hTrigPtEtaPhiPxLFar;
TH3D *hTrigPtEtaPhiPxRClose;
TH3D *hTrigPtEtaPhiPxRFar;

THnSparseD *hDetaDphiVsPtReal;

THnSparseD *hDetaDphiVsPtPxLeftClose;
THnSparseD *hDetaDphiVsPtPxLeftFar;
THnSparseD *hDetaDphiVsPtPxRightClose;
THnSparseD *hDetaDphiVsPtPxRightFar;

THnSparseD *hDetaDphiVsPtMixPxLClose;
THnSparseD *hDetaDphiVsPtMixPxRClose;
THnSparseD *hDetaDphiVsPtMixPxLFar;
THnSparseD *hDetaDphiVsPtMixPxRFar;

const Int_t dimCor = 4;
const Int_t nCorBins[dimCor]={400,400,32,96};//trigPt;AssocPt;dEta;dPhi
const Double_t CorlowBins[dimCor]={0,0,-1.6,-0.5*3.1415926};
const Double_t CorupBins[dimCor]={20,20,1.6,1.5*3.1415926};

//----define-function---
void bookHistograms();
void writeHistograms(char* outFile);
Bool_t passEvent(miniDst* evt, Int_t nTrks);
Bool_t passPion(miniDst* evt,Int_t nTrks);
Bool_t passPionForReal(LorentzVec vec);
Bool_t isPxEvent(Float_t Px,Short_t cen,Float_t PxCut);
Bool_t isPiTrgEvent(miniDst* evt,Int_t nTrks);
void makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
Bool_t GetPtPointer(LorentzVec vecCandidate, Int_t *ptP);
void makeMixPairs( LorentzVec vecTrig, TLorentzVector bufferPi, THnSparseD *hn);
void copyToBufferPxLClose();
void copyToBufferPxLFar();
void copyToBufferPxRClose();
void copyToBufferPxRFar();
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2, Float_t *parL1, Float_t *parL2, Float_t *parR1, Float_t *parR2);

const Bool_t debug = 0;

const float PxLcut[9]={-3.25,-5.75,-9.75,-14.75,-20.75,-26.25,-30.25,-35.75,-37.25};
const float PxRcut[9]={-3.75,-6.25,-10.25,-15.75,-22.25,-28.75,-32.75,-38.75,-40.25};

//---global--var---------------
Float_t  vz = 0.;//tpc vz
Int_t    NTrk = 0;
Short_t  centrality9 = 0.;
Short_t  centrality16 = 0.;
Float_t  PxL = 0.;
Float_t  PxR = 0.;
Bool_t   isPxL;
Bool_t   isPxR;
Bool_t   isPion;
Int_t    iran = 0;
Double_t reWeight = 1.;
//------------------------------

LorentzVec PiCandidate;
FloatVec vecPxL;
FloatVec vecPxR;
LorentzVec PiTrig;
LorentzVec PiAssoc;
//assoc
LorentzVec PxLeftPiClose;
LorentzVec PxRightPiClose;
LorentzVec PxLeftPiFar;
LorentzVec PxRightPiFar;
//trig
LorentzVec PxLTrigPiClose;
LorentzVec PxLTrigPiFar;
LorentzVec PxRTrigPiClose;
LorentzVec PxRTrigPiFar;

//------var-used-for-mixevent-------
const Int_t mCenBins = 9; //default:16
const Int_t mVzBins = 20; //default:20
const Int_t mPtBins = 8; //1-2;2-3;3-4;4-5;5-6,6-7,7-8,>8
const Int_t mMaxEvtsInBuffer = 100; //default:50
const Double_t VzCut = 100.;

Int_t cenPointer, vzPointer;
Int_t ptPointerPxLClose, ptPointerPxLFar, ptPointerPxRClose, ptPointerPxRFar;
Bool_t bufferFullFlag;

Int_t nPxLCloseInBuffer[mCenBins][mVzBins][mPtBins];
Int_t nPxRCloseInBuffer[mCenBins][mVzBins][mPtBins];
Int_t nPxLFarInBuffer[mCenBins][mVzBins][mPtBins];
Int_t nPxRFarInBuffer[mCenBins][mVzBins][mPtBins];

LorentzVec bufferPxLClose[mCenBins][mVzBins][mPtBins][mMaxEvtsInBuffer];
LorentzVec bufferPxRClose[mCenBins][mVzBins][mPtBins][mMaxEvtsInBuffer];
LorentzVec bufferPxLFar[mCenBins][mVzBins][mPtBins][mMaxEvtsInBuffer];
LorentzVec bufferPxRFar[mCenBins][mVzBins][mPtBins][mMaxEvtsInBuffer];
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

    memset(nPxLCloseInBuffer,0,sizeof(nPxLCloseInBuffer));

    miniDst *event = new miniDst(chain);

    Int_t nEvts = chain->GetEntries();
    cout<<nEvts<<" events"<<endl;
    int perevt = nEvts/100;

    //-----timer---------
    TStopwatch timer;
    timer.Start();

    //----loop--events---//
    for(Int_t i=0; i<nEvts; i++) {
        if(i%perevt ==0) cout<<"Looking at "<<((float)i/nEvts)*100<<" % .."<<endl;
        event->GetEntry(i);

        PiTrig.clear();
        PiAssoc.clear();
        PiCandidate.clear();
        vecPxL.clear();
        vecPxR.clear();
        //assoc
        PxLeftPiClose.clear();
        PxLeftPiFar.clear();
        PxRightPiClose.clear();
        PxRightPiFar.clear();
        //trig
        PxLTrigPiClose.clear();
        PxLTrigPiFar.clear();
        PxRTrigPiClose.clear();
        PxRTrigPiFar.clear();

        //--parameters--global-----
        NTrk = event -> mNTrk;
        vz = event -> mTpcVz;
        centrality9 = event -> mCentrality9;
        centrality16 = event -> mCentrality16;

        //-----save--Px--value------------------
        if(!passEvent(event,NTrk)) continue;
        hEvent->Fill(0.5);

        //----2pion--event--cut---------
        Int_t nPi = PiCandidate.size();
        if(debug) cout<<"# of pion candidate: "<<nPi<<endl;
        if(nPi<2) continue;
        hEvent->Fill(1.5);

        //-----choose--px--value---------------
        for(Int_t i=1; i<nPi; i++){
            comparePt(&PiCandidate[0],&PiCandidate[i],&vecPxL[0],&vecPxL[i],&vecPxR[0],&vecPxR[i]);
        }

        //-----get--cen--vz---pointer
        cenPointer = centrality9;
        vzPointer = (Int_t)( (vz+VzCut)/(2*VzCut)*mVzBins );

        isPxL = isPxEvent(vecPxL[0],cenPointer,PxLcut[cenPointer]);
        isPxR = isPxEvent(vecPxR[0],cenPointer,PxRcut[cenPointer]);

        //if(!isPiTrgEvent(event,NTrk)) continue;//Get PiCandidate

        //---Fill event
        if(!isPxL && !isPxR) hEvent->Fill(2.5);
        if(isPxL && !isPxR) hEvent->Fill(3.5);
        if(!isPxL && isPxR) hEvent->Fill(4.5);
        if(isPxL && isPxR) hEvent->Fill(5.5);

        //if(!passPion(event,NTrk)) continue;//Get Pion Candidate

        //*********for check******************
        PiTrig.push_back(PiCandidate[0]);
        for(unsigned int i=1; i<PiCandidate.size(); i++){
            PiAssoc.push_back(PiCandidate[i]);
        }
        makeRealPairs(PiTrig, PiAssoc, hDetaDphiVsPtReal);
        //************************************

        if(!passPionForReal(PiCandidate)) continue;//select by eta range

        //**********PxLClose*****
        unsigned int nTrig = 0, nAssoc = 0; 
        nTrig = PxLTrigPiClose.size();nAssoc = PxLeftPiClose.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(1.5);
            makeRealPairs(PxLTrigPiClose ,PxLeftPiClose,hDetaDphiVsPtPxLeftClose);

            ptPointerPxLClose = -1;
            if(!GetPtPointer(PxLTrigPiClose, &ptPointerPxLClose)) continue;
            if(debug) cout<<"ptPointerPxLClose: "<<ptPointerPxLClose<<endl;
            for(Int_t i=0; i<mPtBins; i++){
                Int_t nEvtInBuf =nPxLCloseInBuffer[cenPointer][vzPointer][ptPointerPxLClose]; 
                for(Int_t j=0; j<nEvtInBuf; j++){
                    Int_t nVecSize = bufferPxLClose[cenPointer][vzPointer][i][j].size();
                    for(Int_t k=0; k<nVecSize; k++){
                        Float_t pt = bufferPxLClose[cenPointer][vzPointer][i][j][k].Pt();
                        Float_t eta = bufferPxLClose[cenPointer][vzPointer][i][j][k].Eta();
                        Float_t phi = bufferPxLClose[cenPointer][vzPointer][i][j][k].Phi();
                        TLorentzVector fourmom(0,0,0,0);
                        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
                        makeMixPairs(PxLTrigPiClose, fourmom, hDetaDphiVsPtMixPxLClose);
                    }
                }
            }
            copyToBufferPxLClose();
            hTrigPtEtaPhiPxLClose->Fill(PxLTrigPiClose[0].Pt(),PxLTrigPiClose[0].Eta(),PxLTrigPiClose[0].Phi());
        }

        //*********PxLFar******
        nTrig = 0; nAssoc = 0;
        nTrig = PxLTrigPiFar.size();nAssoc = PxLeftPiFar.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(2.5);
            makeRealPairs(PxLTrigPiFar ,PxLeftPiFar,hDetaDphiVsPtPxLeftFar);

            ptPointerPxLFar = -1;
            if(!GetPtPointer(PxLTrigPiFar, &ptPointerPxLFar)) continue;
            if(debug) cout<<"ptPointerPxLFar: "<<ptPointerPxLFar<<endl;
            for(Int_t i=0; i<mPtBins; i++){
                Int_t nEvtInBuf =nPxLFarInBuffer[cenPointer][vzPointer][ptPointerPxLFar]; 
                for(Int_t j=0; j<nEvtInBuf; j++){
                    Int_t nVecSize = bufferPxLFar[cenPointer][vzPointer][i][j].size();
                    for(Int_t k=0; k<nVecSize; k++){
                        Float_t pt = bufferPxLFar[cenPointer][vzPointer][i][j][k].Pt();
                        Float_t eta = bufferPxLFar[cenPointer][vzPointer][i][j][k].Eta();
                        Float_t phi = bufferPxLFar[cenPointer][vzPointer][i][j][k].Phi();
                        TLorentzVector fourmom(0,0,0,0);
                        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
                        makeMixPairs(PxLTrigPiFar, fourmom, hDetaDphiVsPtMixPxLFar);
                    }
                }
            }
            copyToBufferPxLFar();
            hTrigPtEtaPhiPxLFar->Fill(PxLTrigPiFar[0].Pt(),PxLTrigPiFar[0].Eta(),PxLTrigPiFar[0].Phi());
        }

        //*********PxRClose******
        nTrig = 0; nAssoc = 0; 
        nTrig = PxRTrigPiClose.size();nAssoc = PxRightPiClose.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(4.5);
            makeRealPairs(PxRTrigPiClose ,PxRightPiClose,hDetaDphiVsPtPxRightClose);

            ptPointerPxRClose = -1;
            if(!GetPtPointer(PxRTrigPiClose, &ptPointerPxRClose)) continue;
            if(debug) cout<<"ptPointerPxRClose: "<<ptPointerPxRClose<<endl;
            for(Int_t i=0; i<mPtBins; i++){
                Int_t nEvtInBuf =nPxRCloseInBuffer[cenPointer][vzPointer][ptPointerPxRClose]; 
                for(Int_t j=0; j<nEvtInBuf; j++){
                    Int_t nVecSize = bufferPxRClose[cenPointer][vzPointer][i][j].size();
                    for(Int_t k=0; k<nVecSize; k++){
                        Float_t pt = bufferPxRClose[cenPointer][vzPointer][i][j][k].Pt();
                        Float_t eta = bufferPxRClose[cenPointer][vzPointer][i][j][k].Eta();
                        Float_t phi = bufferPxRClose[cenPointer][vzPointer][i][j][k].Phi();
                        TLorentzVector fourmom(0,0,0,0);
                        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
                        makeMixPairs(PxRTrigPiClose, fourmom, hDetaDphiVsPtMixPxRClose);
                    }
                }
            }
            copyToBufferPxRClose();
            hTrigPtEtaPhiPxRClose->Fill(PxRTrigPiClose[0].Pt(),PxRTrigPiClose[0].Eta(),PxRTrigPiClose[0].Phi());
        }

        //********PxRFar*******
        nTrig = 0; nAssoc = 0;
        nTrig = PxRTrigPiFar.size();nAssoc = PxRightPiFar.size();
        if( nTrig != 0 && nAssoc != 0 ){
            hnPiTrig->Fill(5.5);
            makeRealPairs(PxRTrigPiFar ,PxRightPiFar,hDetaDphiVsPtPxRightFar);

            ptPointerPxRFar = -1;
            if(!GetPtPointer(PxRTrigPiFar, &ptPointerPxRFar)) continue;
            if(debug) cout<<"ptPointerPxRFar: "<<ptPointerPxRFar<<endl;
            for(Int_t i=0; i<mPtBins; i++){
                Int_t nEvtInBuf =nPxRFarInBuffer[cenPointer][vzPointer][ptPointerPxRFar]; 
                for(Int_t j=0; j<nEvtInBuf; j++){
                    Int_t nVecSize = bufferPxRFar[cenPointer][vzPointer][i][j].size();
                    for(Int_t k=0; k<nVecSize; k++){
                        Float_t pt = bufferPxRFar[cenPointer][vzPointer][i][j][k].Pt();
                        Float_t eta = bufferPxRFar[cenPointer][vzPointer][i][j][k].Eta();
                        Float_t phi = bufferPxRFar[cenPointer][vzPointer][i][j][k].Phi();
                        TLorentzVector fourmom(0,0,0,0);
                        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
                        makeMixPairs(PxRTrigPiFar, fourmom, hDetaDphiVsPtMixPxRFar);
                    }
                }
            }
            copyToBufferPxRFar();
            hTrigPtEtaPhiPxRFar->Fill(PxRTrigPiFar[0].Pt(),PxRTrigPiFar[0].Eta(),PxRTrigPiFar[0].Phi());
        }

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
//------------------------------------------
Bool_t passEvent(miniDst* evt, Int_t nTrks)
{
    //--for Check----
    Double_t vpdVz = evt->mVpdVz;
    Double_t tpcVx = evt->mTpcVx;
    Double_t tpcVy = evt->mTpcVy;
    Double_t tpcVz = evt->mTpcVz;
    Double_t vzDiff = tpcVz - vpdVz;
    Double_t grefMult = evt->mGRefMult;
    Double_t grefMultCorr = evt->mGRefMultCorr;

    hVertexYvsX->Fill(tpcVx,tpcVy);
    hVpdVzvsTpcVz->Fill(vpdVz,tpcVz);
    hVzDiff->Fill(vzDiff);
    hRefMult->Fill(grefMult);
    hGRefMultvsGRefMultCorr->Fill(grefMultCorr,grefMult);

    //--pass PionCandidate and px value
    for(Int_t i=0;i<nTrks;i++){
        Float_t pt = evt->mTrkPt[i];
        Float_t eta = evt->mTrkEta[i];
        Float_t phi = evt->mTrkPhi[i];
        Float_t sigmaPi = evt->mnSigmaPi[i];
        if(!(sigmaPi>0. && sigmaPi<3.)) continue;//choose pion and save
        Float_t pxL = evt -> mPxL[i];
        vecPxL.push_back(pxL);
        Float_t pxR = evt -> mPxR[i];
        vecPxR.push_back(pxR);
        TLorentzVector fourmom(0,0,0,0);
        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
        PiCandidate.push_back(fourmom);
    }
    return kTRUE;
}
//-------------------------------------------
Bool_t passPion(miniDst* evt, Int_t nTrks )
{
    for(Int_t i=0;i<nTrks;i++){
        Float_t pt = evt->mTrkPt[i];//trigger particle
        Float_t eta = evt->mTrkEta[i];
        Float_t phi = evt->mTrkPhi[i];
        Float_t sigmaPi = evt->mnSigmaPi[i];
        if(!(sigmaPi>0. && sigmaPi<3.)) continue;
        TLorentzVector fourmom(0,0,0,0);
        fourmom.SetPtEtaPhiM(pt,eta,phi,Mpi);
        PiCandidate.push_back(fourmom);
    }
    if(PiCandidate.size()>1) return kTRUE;
    else return kFALSE;
}
//-------------------------------------------
Bool_t passPionForReal(LorentzVec vec)
{
    Int_t nvec = vec.size();
    if(nvec<2) {
        if(debug) cout<<"PiCandidate error nPiCandidate: "<<nvec<<endl;
        return kFALSE;
    }
    if(debug) cout << "# of pion candidate: "<< nvec << endl;

    for(Int_t i=1; i<nvec; i++){
        comparePt(&vec[0],&vec[i]);
    }
    if(debug) cout<<"PassPionForReal: trigPt: "<<vec[0].Pt()<<endl;

    Int_t nPxLClose = 0;
    Int_t nPxLFar = 0;
    if(isPxL){
        hEvent->Fill(6.5);
        hnPiTrig->Fill(0.5);
        for(Int_t i=1; i<nvec; i++){
            Float_t eta = vec[i].Eta();
            if(eta<0 && eta>-0.5) {
                PxLeftPiClose.push_back(vec[i]);
                nPxLClose++;
            }
            if(eta>0 && eta<0.5) {
                PxLeftPiFar.push_back(vec[i]);
                nPxLFar++;
            }
        }
        if( nPxLClose>0 && nPxLClose==(Int_t)PxLeftPiClose.size() ) 
        {PxLTrigPiClose.push_back(vec[0]);}
        if( nPxLFar>0 && nPxLFar==(Int_t)PxLeftPiFar.size() ) 
            PxLTrigPiFar.push_back(vec[0]);
    }

    Int_t nPxRClose = 0;
    Int_t nPxRFar = 0;
    if(isPxR){
        hEvent->Fill(7.5);
        hnPiTrig->Fill(3.5);
        for(Int_t i=1;i<nvec ; i++){
            Float_t eta = vec[i].Eta();
            if(eta<0 && eta>-0.5){
                PxRightPiFar.push_back(vec[i]);
                nPxRFar++;
            }
            if(debug && nPxRFar==0) cout<<"PxR not far eta: "<<eta<<endl;
            if(eta>0 && eta<0.5){
                PxRightPiClose.push_back(vec[i]);
                nPxRClose++;
            }
            if(debug && nPxRClose==0 && nPxRFar==0 ) cout<<"PxR not both close and far eta: "<<eta<<endl;
        }
        if( nPxRClose>0 && nPxRClose==(Int_t)PxRightPiClose.size() ) 
        {PxRTrigPiClose.push_back(vec[0]);}
        if( nPxRFar>0 && nPxRFar==(Int_t)PxRightPiFar.size() ) 
        {PxRTrigPiFar.push_back(vec[0]);}
    }
    return kTRUE;
}
//-------------------------------------------
Bool_t isPxEvent(Float_t Px,Short_t cen,Float_t PxCut)
{
    if(Px < PxCut) {
        if(debug) cout<<"Px: "<<Px<<" cen: "<<cen<<" PxCut: "<<PxCut<<endl;
        return kTRUE;
    }
    else return kFALSE;
}
//-------------------------------------------
Bool_t is2PiEvent( miniDst* evt, Int_t nTrks )
{
    Int_t nPi = 0;
    for(Int_t i=0;i<nTrks;i++){
        Float_t sigmaPi = evt->mnSigmaPi[i];
        if(!(sigmaPi>0. && sigmaPi<3.)) continue;
        nPi++;
    }
    if(nPi>1) return kTRUE;
    else return kFALSE;
}
//-------------------------------------------
Bool_t isPiTrgEvent( miniDst* evt, Int_t nTrks )
{
    if(!is2PiEvent(evt,nTrks)) return kFALSE;
    Float_t pt[100] = {0.};
    Float_t sigPi[100] = {0.};
    for(Int_t i=0;i<nTrks;i++){
        pt[i] = evt->mTrkPt[i];  
        sigPi[i] = evt->mnSigmaPi[i];  
    }
    for(Int_t i=1;i<nTrks;i++){
        Float_t ptp = 0.,sigPip = 0.; 
        if(pt[0]<pt[i]){
            ptp = pt[0];pt[0]=pt[i];pt[i]=ptp;
            sigPip = sigPi[0];sigPi[0]=sigPi[i];sigPi[i]=sigPip;
        }
    }
    if(sigPi[0]>0. && sigPi[0]<3.) return kTRUE;
    else return kFALSE;
}
//-------------------------------------------
void makeRealPairs(LorentzVec vecTrig, LorentzVec vecAssoc,THnSparseD *hn)
{
    Float_t TrigPt = vecTrig[0].Pt();
    Float_t TrigEta = vecTrig[0].Eta();
    Float_t TrigPhi = vecTrig[0].Phi();

    Int_t nVec = vecAssoc.size();
    for(Int_t i=0; i<nVec; i++){
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
}
//-------------------------------------------
Bool_t GetPtPointer(LorentzVec vecCandidate, Int_t *ptPointer)
{//to select trigger particle 
    Int_t nVec = vecCandidate.size();
    if(nVec<1)     return kFALSE;
    if(nVec!=1){
        for(Int_t i=0;i<nVec;i++){
            comparePt(&vecCandidate[0],&vecCandidate[i]);
        }}
    Float_t candidatePt = vecCandidate[0].Pt();

    Int_t ptP = -1;
    if(1.<candidatePt && candidatePt<2.) ptP=0;
    else if(2.<candidatePt && candidatePt<3.) ptP=1;
    else if(3.<candidatePt && candidatePt<4.) ptP=2;
    else if(4.<candidatePt && candidatePt<5.) ptP=3;
    else if(5.<candidatePt && candidatePt<6.) ptP=4;
    else if(6.<candidatePt && candidatePt<7.) ptP=5;
    else if(7.<candidatePt && candidatePt<8.) ptP=6;
    else if(8.<candidatePt) ptP=7;

    *ptPointer = ptP;
    return kTRUE;
} 
//-------------------------------------------
void makeMixPairs( LorentzVec vecTrig, TLorentzVector bufferPi, THnSparseD *hn)
{
    Float_t TrigPt = vecTrig[0].Pt();
    Float_t TrigEta = vecTrig[0].Eta();
    Float_t TrigPhi = vecTrig[0].Phi();

    Float_t BufferPt = bufferPi.Pt();
    Float_t BufferEta = bufferPi.Eta();
    Float_t BufferPhi = bufferPi.Phi();

    Float_t dEta = BufferEta - TrigEta;
    Float_t dPhi = BufferPhi - TrigPhi;
    if(dPhi < -TMath::Pi()/2.) dPhi += 2*TMath::Pi();
    if(dPhi > 3*TMath::Pi()/2.) dPhi -= 2*TMath::Pi();
    Double_t fill[]={TrigPt,BufferPt,dEta,dPhi};
    hn->Fill(fill);
}
//-------------------------------------------
void copyToBufferPxLClose()
{
    if(nPxLCloseInBuffer[cenPointer][vzPointer][ptPointerPxLClose]>=mMaxEvtsInBuffer){
        bufferFullFlag = kTRUE;
    }//full flag
    TRandom3 *gRandom = new TRandom3(iran++);
    Int_t evtPointer = -1;
    if( bufferFullFlag ) evtPointer = (Int_t)gRandom->Uniform(0,mMaxEvtsInBuffer-1.e-6);
    else evtPointer = nPxLCloseInBuffer[cenPointer][vzPointer][ptPointerPxLClose];
    delete gRandom;

    bufferPxLClose[cenPointer][vzPointer][ptPointerPxLClose][evtPointer].clear();

    bufferPxLClose[cenPointer][vzPointer][ptPointerPxLClose][evtPointer].push_back(PxLTrigPiClose[0]);
    Int_t nVec = PxLeftPiClose.size();
    for(Int_t i=0;i<nVec;i++){
        bufferPxLClose[cenPointer][vzPointer][ptPointerPxLClose][evtPointer].push_back(PxLeftPiClose[i]);
    }
    if(nPxLCloseInBuffer[cenPointer][vzPointer][ptPointerPxLClose] < mMaxEvtsInBuffer){
        nPxLCloseInBuffer[cenPointer][vzPointer][ptPointerPxLClose]++;
    }
}
//-------------------------------------------
void copyToBufferPxLFar()
{
    if(nPxLFarInBuffer[cenPointer][vzPointer][ptPointerPxLFar]>=mMaxEvtsInBuffer){
        bufferFullFlag = kTRUE;
    }//full flag
    TRandom3 *gRandom = new TRandom3(iran++);
    Int_t evtPointer = -1;
    if( bufferFullFlag ) evtPointer = (Int_t)gRandom->Uniform(0,mMaxEvtsInBuffer-1.e-6);
    else evtPointer = nPxLFarInBuffer[cenPointer][vzPointer][ptPointerPxLFar];
    delete gRandom;

    bufferPxLFar[cenPointer][vzPointer][ptPointerPxLFar][evtPointer].clear();

    bufferPxLFar[cenPointer][vzPointer][ptPointerPxLFar][evtPointer].push_back(PxLTrigPiFar[0]);
    Int_t nVec = PxLeftPiFar.size();
    for(Int_t i=0;i<nVec;i++){
        bufferPxLFar[cenPointer][vzPointer][ptPointerPxLFar][evtPointer].push_back(PxLeftPiFar[i]);
    }
    if(nPxLFarInBuffer[cenPointer][vzPointer][ptPointerPxLFar] < mMaxEvtsInBuffer){
        nPxLFarInBuffer[cenPointer][vzPointer][ptPointerPxLFar]++;
    }
}
//-------------------------------------------
void copyToBufferPxRClose()
{
    if(nPxRCloseInBuffer[cenPointer][vzPointer][ptPointerPxRClose]>=mMaxEvtsInBuffer){
        bufferFullFlag = kTRUE;
    }//full flag
    TRandom3 *gRandom = new TRandom3(iran++);
    Int_t evtPointer = -1;
    if( bufferFullFlag ) evtPointer = (Int_t)gRandom->Uniform(0,mMaxEvtsInBuffer-1.e-6);
    else evtPointer = nPxRCloseInBuffer[cenPointer][vzPointer][ptPointerPxRClose];
    delete gRandom;

    bufferPxRClose[cenPointer][vzPointer][ptPointerPxRClose][evtPointer].clear();

    bufferPxRClose[cenPointer][vzPointer][ptPointerPxRClose][evtPointer].push_back(PxRTrigPiClose[0]);
    Int_t nVec = PxRightPiClose.size();
    for(Int_t i=0;i<nVec;i++){
        bufferPxRClose[cenPointer][vzPointer][ptPointerPxRClose][evtPointer].push_back(PxRightPiClose[i]);
    }
    if(nPxRCloseInBuffer[cenPointer][vzPointer][ptPointerPxRClose] < mMaxEvtsInBuffer){
        nPxRCloseInBuffer[cenPointer][vzPointer][ptPointerPxRClose]++;
    }
}
//-------------------------------------------
void copyToBufferPxRFar()
{
    if(nPxRFarInBuffer[cenPointer][vzPointer][ptPointerPxRFar]>=mMaxEvtsInBuffer){
        bufferFullFlag = kTRUE;
    }//full flag
    TRandom3 *gRandom = new TRandom3(iran++);
    Int_t evtPointer = -1;
    if( bufferFullFlag ) evtPointer = (Int_t)gRandom->Uniform(0,mMaxEvtsInBuffer-1.e-6);
    else evtPointer = nPxRFarInBuffer[cenPointer][vzPointer][ptPointerPxRFar];
    delete gRandom;

    bufferPxRFar[cenPointer][vzPointer][ptPointerPxRFar][evtPointer].clear();

    bufferPxRFar[cenPointer][vzPointer][ptPointerPxRFar][evtPointer].push_back(PxRTrigPiFar[0]);
    Int_t nVec = PxRightPiFar.size();
    for(Int_t i=0;i<nVec;i++){
        bufferPxRFar[cenPointer][vzPointer][ptPointerPxRFar][evtPointer].push_back(PxRightPiFar[i]);
    }
    if(nPxRFarInBuffer[cenPointer][vzPointer][ptPointerPxRFar] < mMaxEvtsInBuffer){
        nPxRFarInBuffer[cenPointer][vzPointer][ptPointerPxRFar]++;
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
    hEvent = new TH1D("hEvent","# of event",8,0,8);
    hEvent->GetXaxis()->SetBinLabel(1,"all dimuon event");
    hEvent->GetXaxis()->SetBinLabel(2,"2Pion cut");
    hEvent->GetXaxis()->SetBinLabel(3,"reject event");
    hEvent->GetXaxis()->SetBinLabel(4,"PxL only");
    hEvent->GetXaxis()->SetBinLabel(5,"PxR only");
    hEvent->GetXaxis()->SetBinLabel(6,"PxL and PxR");
    hEvent->GetXaxis()->SetBinLabel(7,"PxL");
    hEvent->GetXaxis()->SetBinLabel(8,"PxR");

    hnPiTrig = new TH1D("hnPiTrig","# of Triggered pion",10,0,10);
    hnPiTrig->GetXaxis()->SetBinLabel(1,"nPxLeftClose");
    hnPiTrig->GetXaxis()->SetBinLabel(2,"nPxLeftClose");
    hnPiTrig->GetXaxis()->SetBinLabel(3,"nPxLeftFar");
    hnPiTrig->GetXaxis()->SetBinLabel(4,"nPxRightClose");
    hnPiTrig->GetXaxis()->SetBinLabel(5,"nPxRightClose");
    hnPiTrig->GetXaxis()->SetBinLabel(6,"nPxRightFar");
    hnPiTrig->GetXaxis()->SetBinLabel(7,"nBothLeftRight");
    hnPiTrig->GetXaxis()->SetBinLabel(8,"nEvents");

    hTrigPtEtaPhiPxLClose = new TH3D("hTrigPtEtaPhiPxLClose","PxL close-region triggered ;p_{T trig} (GeV/c);#eta;#phi",400,0,20,20,-1.,1.,96,-3.1415926,3.1415926);
    hTrigPtEtaPhiPxLFar = new TH3D("hTrigPtEtaPhiPxLFar","PxL close-region triggered ;p_{T trig} (GeV/c);#eta;#phi",400,0,20,20,-1.,1.,96,-3.1415926,3.1415926);
    hTrigPtEtaPhiPxRClose = new TH3D("hTrigPtEtaPhiPxRClose","PxR close-region triggered ;p_{T trig} (GeV/c);#eta;#phi",400,0,20,20,-1.,1.,96,-3.1415926,3.1415926);
    hTrigPtEtaPhiPxRFar = new TH3D("hTrigPtEtaPhiPxRFar","PxR close-region triggered ;p_{T trig} (GeV/c);#eta;#phi",400,0,20,20,-1.,1.,96,-3.1415926,3.1415926);

    hDetaDphiVsPtReal = new THnSparseD("hDetaDphiVsPtReal","for Check & Compare ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);

    hDetaDphiVsPtPxLeftClose = new THnSparseD("hDetaDphiVsPtPxLeftClose","Px  -1.0 < #eta < -0.5  close region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxLeftFar = new THnSparseD("hDetaDphiVsPtPxLeftFar"," Px  -1.0 < #eta < -0.5 far region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxRightClose = new THnSparseD("hDetaDphiVsPtPxRightClose"," Px  0.5 < #eta < 1.0  close region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxRightFar = new THnSparseD("hDetaDphiVsPtPxRightFar"," Px  0.5 < #eta < 1.0  far region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);

    hDetaDphiVsPtMixPxLClose = new THnSparseD("hDetaDphiVsPtMixPxLClose","mixevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtMixPxRClose = new THnSparseD("hDetaDphiVsPtMixPxRClose","mixevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtMixPxLFar = new THnSparseD("hDetaDphiVsPtMixPxLFar","mixevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtMixPxRFar = new THnSparseD("hDetaDphiVsPtMixPxRFar","mixevent ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
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

    hEvent->Write();
    hnPiTrig->Write();

hTrigPtEtaPhiPxLClose->Write();
hTrigPtEtaPhiPxLFar->Write();
hTrigPtEtaPhiPxRClose->Write();
hTrigPtEtaPhiPxRFar->Write();

    hDetaDphiVsPtReal -> Write();

    hDetaDphiVsPtPxLeftClose -> Write();
    hDetaDphiVsPtPxLeftFar -> Write();
    hDetaDphiVsPtPxRightClose -> Write();
    hDetaDphiVsPtPxRightFar -> Write();

    hDetaDphiVsPtMixPxLClose -> Write();
    hDetaDphiVsPtMixPxRClose -> Write();
    hDetaDphiVsPtMixPxLFar -> Write();
    hDetaDphiVsPtMixPxRFar -> Write();
}
