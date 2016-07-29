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
TH1D *hTrigPt;
TH1D *hLeftTrigPt;
TH1D *hRightTrigPt;
TH1D *hAssocPt;
TH1D *hLeftAssocPt;
TH1D *hRightAssocPt;

THnSparseD *hDetaDphiVsPtPxLeftClose;
THnSparseD *hDetaDphiVsPtPxLeftFar;
THnSparseD *hDetaDphiVsPtPxRightClose;
THnSparseD *hDetaDphiVsPtPxRightFar;

THnSparseD *hDetaDphiVsPtMix;

const Int_t dimCor = 4;
const Int_t nCorBins[dimCor]={400,400,32,96};//trigPt;AssocPt;dEta;dPhi
const Double_t CorlowBins[dimCor]={0,0,-1.6,-0.5*3.1415926};
const Double_t CorupBins[dimCor]={20,20,1.6,1.5*3.1415926};

//----define-function---
void bookHistograms();
void writeHistograms(char* outFile);
Bool_t passEvent(miniDst* evt);
Bool_t passPion(miniDst* evt,Int_t nTrks);
Bool_t passPionForReal(LorentzVec vec);
Bool_t isPxEvent(Float_t Px,Short_t cen,Float_t PxCut);
Bool_t isPiTrgEvent(miniDst* evt,Int_t nTrks);
Bool_t GetPtPointer(LorentzVec vecCandidate);
Bool_t makeMixPairs(LorentzVec vecTrig,THnSparseD *hn);
void copyToBuffer(LorentzVec vecCandidate);
void makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void copyToBuffer(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2);

const Bool_t debug = 0;

const float PxLcut[9]={-3.25,-5.75,-9.75,-14.75,-20.75,-27.75,-31.75,-37.25,-38.25};
const float PxRcut[9]={-3.75,-6.25,-10.25,-15.75,-22.25,-29.75,-34.25,-40.25,-41.75};

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
LorentzVec PiTrig;
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
    for(Int_t i=0; i<nEvts; i++) {
        if(i%perevt ==0) cout<<"Looking at "<<((float)i/nEvts)*100<<" % .."<<endl;
        event->GetEntry(i);

        PiTrig.clear();
        PiCandidate.clear();

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
        PxL = event -> mPxL;
        PxR = event -> mPxR;
        if(debug)cout<<"NTrk: "<<NTrk;
        if(debug)cout<<"  TpcVz: "<<vz;
        if(debug)cout<<"  centrality9: "<<centrality9;
        if(debug)cout<<"  centrality16: "<<centrality16;
        if(debug)cout<<"  PxL: "<<PxL;
        if(debug)cout<<"  PxR: "<<PxR<<endl;
        //----------------------------
        if(!passEvent(event)) continue;
        hEvent->Fill(0.5);

        //get cen vz  pointer
        cenPointer = centrality9;
        vzPointer = (Int_t)( (vz+VzCut)/(2*VzCut)*mVzBins );
        if(debug)cout<<"vzPointer:  "<<vzPointer<<endl;

        isPxL = isPxEvent(PxL,cenPointer,PxLcut[cenPointer]);
        isPxR = isPxEvent(PxR,cenPointer,PxRcut[cenPointer]);

        if(!isPiTrgEvent(event,NTrk)) continue;//Get PiCandidate
        hEvent->Fill(1.5);
        if(!isPxL && !isPxR) hEvent->Fill(2.5);
        if(isPxL && !isPxR) hEvent->Fill(3.5);
        if(!isPxL && isPxR) hEvent->Fill(4.5);
        if(isPxL && isPxR) hEvent->Fill(5.5);

        if(!passPion(event,NTrk)) continue;//Get Pion Candidate
        if(!passPionForReal(PiCandidate)) continue;//select by eta range

        unsigned int nTrig = 0, nAssoc = 0; 
        nTrig = PxLTrigPiClose.size();nAssoc = PxLeftPiClose.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(1.5);
            makeRealPairs(PxLTrigPiClose ,PxLeftPiClose,hDetaDphiVsPtPxLeftClose);}

        nTrig = 0; nAssoc = 0;
        nTrig = PxLTrigPiFar.size();nAssoc = PxLeftPiFar.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(2.5);
            makeRealPairs(PxLTrigPiFar ,PxLeftPiFar,hDetaDphiVsPtPxLeftFar);}

        nTrig = 0; nAssoc = 0; 
        nTrig = PxRTrigPiClose.size();nAssoc = PxRightPiClose.size();
        if( nTrig != 0 && nAssoc != 0 ){ 
            hnPiTrig->Fill(4.5);
            makeRealPairs(PxRTrigPiClose ,PxRightPiClose,hDetaDphiVsPtPxRightClose);}

        nTrig = 0; nAssoc = 0;
        nTrig = PxRTrigPiFar.size();nAssoc = PxRightPiFar.size();
        if( nTrig != 0 && nAssoc != 0 ){
            hnPiTrig->Fill(5.5);
            makeRealPairs(PxRTrigPiFar ,PxRightPiFar,hDetaDphiVsPtPxRightFar);}

        //GetPtPointer(PiCandidate);
        //if(debug) cout<<"ptPointer: "<<ptPointer<<endl;


        //hTrigPt->Fill(PiCandidate[0].Pt());
        //for(unsigned int i=1;i<PiCandidate.size();i++) hAssocPt->Fill(PiCandidate[i].Pt());//associate
        //hnPiTrig->Fill(7.5);//trigger

        //if(debug) cout<<"# of Trigger Particle "<<PiTrig.size()<<endl;
        //if(isPxL){
        //    hLeftTrigPt->Fill(PiCandidate[0].Pt());
        //    hnPiTrig->Fill(0.5);
        //    if(makeRealPairs(PiTrig,PxLeftPiClose,hDetaDphiVsPtPxLeftClose)) 
        //        hnPiTrig->Fill(1.5);
        //    if(PxLeftPiClose.size()!=0){ for(unsigned int i=0;i<PxLeftPiClose.size();i++) hLeftAssocPt->Fill(PxLeftPiClose[i].Pt());}
        //    if(makeRealPairs(PiTrig,PxLeftPiFar,hDetaDphiVsPtPxLeftFar)) 
        //        hnPiTrig->Fill(2.5);
        //    if(PxLeftPiFar.size()!=0){ for(unsigned int i=0;i<PxLeftPiFar.size();i++) hLeftAssocPt->Fill(PxLeftPiFar[i].Pt());}
        //}
        //if(isPxR){
        //    hRightTrigPt->Fill(PiCandidate[0].Pt());
        //    hnPiTrig->Fill(3.5);
        //    if(makeRealPairs(PiTrig,PxRightPiClose,hDetaDphiVsPtPxRightClose)) 
        //        hnPiTrig->Fill(4.5);
        //    if(PxRightPiClose.size()!=0){ for(unsigned int i=0;i<PxRightPiClose.size();i++) hRightAssocPt->Fill(PxRightPiClose[i].Pt());}
        //    if(makeRealPairs(PiTrig,PxRightPiFar,hDetaDphiVsPtPxRightFar)) 
        //        hnPiTrig->Fill(5.5);
        //    if(PxRightPiFar.size()!=0){ for(unsigned int i=0;i<PxRightPiFar.size();i++) hRightAssocPt->Fill(PxRightPiFar[i].Pt());}
        //}
        //if(isPxL && isPxR){
        //    hnPiTrig->Fill(6.5);
        //}

        //if(debug)cout<<"mix event"<<endl;
        //makeMixPairs(PiTrig,hDetaDphiVsPtMix);
        //copyToBuffer(PiCandidate);

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

    Int_t nPxLClose = 0;
    Int_t nPxLFar = 0;
    if(isPxL){
        hEvent->Fill(6.5);
        for(Int_t i=1; i<nvec; i++){
            Float_t eta = vec[i].Eta();
            if(debug) cout<<eta<<endl;
            if(eta<0 && eta>-0.5) {
                if(debug) cout<<"isPxLClose eta: "<<eta;
                PxLeftPiClose.push_back(vec[i]);
                nPxLClose++;
            }
            else if(eta>0 && eta<0.5) {
                if(debug) cout<<"isPxLFar eta: "<<eta;
                PxLeftPiFar.push_back(vec[i]);
                nPxLFar++;
            }
        }
        if( nPxLClose>0 && nPxLClose==(Int_t)PxLeftPiClose.size() ) 
        {PxLTrigPiClose.push_back(vec[0]);hnPiTrig->Fill(0.5);}
        if(nPxLFar>0 && nPxLFar==(Int_t)PxLeftPiFar.size() ) 
            PxLTrigPiFar.push_back(vec[0]);
        if(debug){
            cout<<"passPionForReal isPxL: "<<isPxL
                <<" ; nClose: "<<nPxLClose<<" ; "<<PxLeftPiClose.size()
                <<" ; TrigClose: "<<PxLeftPiClose.size()
                <<" ; nFar: "<<nPxLFar<<" ; "<<PxLeftPiFar.size()
                <<" ; TrigFar: "<<PxLeftPiFar.size()
                <<" ; "<<endl;
        }
    }

    Int_t nPxRClose = 0;
    Int_t nPxRFar = 0;
    if(isPxR){
        hEvent->Fill(7.5);
        for(Int_t i=1;i<nvec ; i++){
            Float_t eta = vec[i].Eta();
            if(debug) cout<<eta<<endl;
            if(eta<0 && eta>-0.5){
                if(debug) cout<<"isPxRFar eta: "<<eta;
                PxRightPiFar.push_back(vec[i]);
                nPxRClose++;
            }
            else if(eta>0 && eta<0.5){
                if(debug) cout<<"isPxRClose eta: "<<eta;
                PxRightPiClose.push_back(vec[i]);
                nPxRFar++;
            }
        }
        if( nPxRClose>0 && nPxRClose==(Int_t)PxRightPiClose.size() ) 
        {PxRTrigPiClose.push_back(vec[0]);hnPiTrig->Fill(3.5);}
        if(nPxRFar>0 && nPxRFar==(Int_t)PxRightPiFar.size() ) 
            PxRTrigPiFar.push_back(vec[0]);
        if(debug){
            cout<<"passPionForReal isPxR: "<<isPxR
                <<" ; nClose: "<<nPxRClose<<" ; "<<PxRightPiClose.size()
                <<" ; TrigClose: "<<PxRightPiClose.size()
                <<" ; nFar: "<<nPxRFar<<" ; "<<PxRightPiFar.size()
                <<" ; TrigFar: "<<PxRightPiFar.size()
                <<" ; "<<endl;
        }
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
void makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn)
{
    unsigned int nVecAssoc = vecAssoc.size();
    Float_t TrigPt = vecTrig[0].Pt();
    Float_t TrigEta = vecTrig[0].Eta();
    Float_t TrigPhi = vecTrig[0].Phi();

    for(unsigned int i =0;i<nVecAssoc;i++){
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
    hnPiTrig->GetXaxis()->SetBinLabel(1,"nPxLeft");
    hnPiTrig->GetXaxis()->SetBinLabel(2,"nPxLeftClose");
    hnPiTrig->GetXaxis()->SetBinLabel(3,"nPxLeftFar");
    hnPiTrig->GetXaxis()->SetBinLabel(4,"nPxRight");
    hnPiTrig->GetXaxis()->SetBinLabel(5,"nPxRightClose");
    hnPiTrig->GetXaxis()->SetBinLabel(6,"nPxRightFar");
    hnPiTrig->GetXaxis()->SetBinLabel(7,"nBothLeftRight");
    hnPiTrig->GetXaxis()->SetBinLabel(8,"nEvents");

    hTrigPt = new TH1D("hTrigPt"," p_{T}^{trig} event level ;p_{T}^{trig} (GeV/c)",400,0,20);
    hLeftTrigPt = new TH1D("hLeftTrigPt"," -1.0 < #eta < -0.5  p_{T}^{trig} ;p_{T}^{trig} (GeV/c)",400,0,20);
    hRightTrigPt = new TH1D("hRightTrigPt"," 0.5 < #eta < 1.0  p_{T}^{trig} ;p_{T}^{trig} (GeV/c)",400,0,20);
    hAssocPt = new TH1D("hAssocPt","hAssocPt;p_{T}^{assoc} (GeV/c)",400,0,20);
    hLeftAssocPt = new TH1D("hLeftAssocPt"," -1.0 < #eta < -0.5  p_{T}^{assoc} ;p_{T}^{assoc} (GeV/c)",400,0,20);
    hRightAssocPt = new TH1D("hRightAssocPt"," 0.5 < #eta < 1.0  p_{T}^{assoc} ;p_{T}^{assoc} (GeV/c)",400,0,20);

    hDetaDphiVsPtPxLeftClose = new THnSparseD("hDetaDphiVsPtPxLeftClose","Px  -1.0 < #eta < -0.5  close region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxLeftFar = new THnSparseD("hDetaDphiVsPtPxLeftFar"," Px  -1.0 < #eta < -0.5 far region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxRightClose = new THnSparseD("hDetaDphiVsPtPxRightClose"," Px  0.5 < #eta < 1.0  close region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);
    hDetaDphiVsPtPxRightFar = new THnSparseD("hDetaDphiVsPtPxRightFar"," Px  0.5 < #eta < 1.0  far region ;p_{T trig} (GeV/c);p_{T assoc} (GeV/c);#Delta#eta;#Delta#phi",dimCor,nCorBins,CorlowBins,CorupBins);

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
    hEvent->Write();

    hnPiTrig->Write();
    hTrigPt->Write();
    hLeftTrigPt->Write();
    hRightTrigPt->Write();
    hAssocPt->Write();
    hLeftAssocPt->Write();
    hRightAssocPt->Write();

    hDetaDphiVsPtPxLeftClose -> Write();
    hDetaDphiVsPtPxLeftFar -> Write();
    hDetaDphiVsPtPxRightClose -> Write();
    hDetaDphiVsPtPxRightFar -> Write();

    hDetaDphiVsPtMix -> Write();
}
