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
TH1D *hCentrality;

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
Bool_t passPion(miniDst* evt);
Bool_t passPionForReal(LorentzVec vec);
Int_t GetCenBinPointer(Short_t cen);
Bool_t IsPxEvent(Float_t Px,Short_t cen,Float_t PxCut);
Bool_t GetPtPointer(LorentzVec vecCandidate);
Bool_t makeMixPairs(LorentzVec vecTrig,THnSparseD *hn);
void copyToBuffer(LorentzVec vecCandidate);
Bool_t makeRealPairs(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void copyToBuffer(LorentzVec vecTrig,LorentzVec vecAssoc,THnSparseD *hn);
void comparePt(TLorentzVector* vec1 ,TLorentzVector* vec2);

const Bool_t debug = 1;

const float PxLcut[9]={-3.25,-5.75,-9.75,-14.75,-20.75,-27.75,-31.75,-37.25,-38.25};
const float PxRcut[9]={-3.75,-6.25,-10.25,-15.75,-22.25,-29.75,-34.25,-40.25,-41.75};

//---global--var---------------
Float_t  vz = 0.;//tpc vz
Int_t    NPiTrks = 0;
Short_t  centrality = 0.;
Float_t  PxL = 0.;
Float_t  PxR = 0.;
Bool_t   IsPxL;
Bool_t   IsPxR;
Int_t    iran = 0;
Double_t reWeight = 1.;
//------------------------------

LorentzVec PiCandidate;
LorentzVec PiTrig;
LorentzVec PxLeftPiClose;
LorentzVec PxRightPiClose;
LorentzVec PxLeftPiFar;
LorentzVec PxRightPiFar;

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
        if(debug)cout<<i<<endl;

        PiTrig.clear();
        PiCandidate.clear();
        PxLeftPiClose.clear();
        PxLeftPiFar.clear();
        PxRightPiClose.clear();
        PxRightPiFar.clear();

        //--parameters--global-----
        NPiTrks = event -> mNPiTrk;
        vz = event -> mTpcVz;
        centrality = event -> mCentrality;
        PxL = event -> mPxL;
        PxR = event -> mPxR;
        if(debug)cout<<"NPiTrks: "<<NPiTrks;
        if(debug)cout<<"  TpcVz: "<<vz;
        if(debug)cout<<"  centrality: "<<centrality;
        if(debug)cout<<"  PxL: "<<PxL;
        if(debug)cout<<"  PxR: "<<PxR<<endl;
        //----------------------------

        if(!passEvent(event)) continue;
        Int_t cenbin = GetCenBinPointer(centrality);
        //get cen vz  pointer
        cenPointer = cenbin;
        vzPointer = (Int_t)( (vz+VzCut)/(2*VzCut)*mVzBins );
        if(debug)cout<<"centrality: "<<centrality<<"  cenpointer: "<<cenPointer<<"  vzPointer:  "<<vzPointer<<endl;

        IsPxL = IsPxEvent(PxL,centrality,PxLcut[cenbin]);
        if(debug && IsPxL) cout<<"IsPxL event"<<endl;
        IsPxR = IsPxEvent(PxR,centrality,PxRcut[cenbin]);
        if(debug && IsPxR) cout<<"IsPxR event"<<endl;

        if(!passPion(event)) continue;//Get PiCandidate

        GetPtPointer(PiCandidate);
        if(debug) cout<<"ptPointer: "<<ptPointer<<endl;

        if(!passPionForReal(PiCandidate)) continue;//Get Trigger&Assoc

        hTrigPt->Fill(PiCandidate[0].Pt());
        for(unsigned int i=1;i<PiCandidate.size();i++) hAssocPt->Fill(PiCandidate[i].Pt());//associate
        hnPiTrig->Fill(7.5);//trigger

        if(debug) cout<<"# of Trigger Particle "<<PiTrig.size()<<endl;
        if(IsPxL){
            hLeftTrigPt->Fill(PiCandidate[0].Pt());
            hnPiTrig->Fill(0.5);
            if(makeRealPairs(PiTrig,PxLeftPiClose,hDetaDphiVsPtPxLeftClose)) 
                hnPiTrig->Fill(1.5);
            if(PxLeftPiClose.size()!=0){ for(unsigned int i=0;i<PxLeftPiClose.size();i++) hLeftAssocPt->Fill(PxLeftPiClose[i].Pt());}
            if(makeRealPairs(PiTrig,PxLeftPiFar,hDetaDphiVsPtPxLeftFar)) 
                hnPiTrig->Fill(2.5);
            if(PxLeftPiFar.size()!=0){ for(unsigned int i=0;i<PxLeftPiFar.size();i++) hLeftAssocPt->Fill(PxLeftPiFar[i].Pt());}
        }
        if(IsPxR){
            hRightTrigPt->Fill(PiCandidate[0].Pt());
            hnPiTrig->Fill(3.5);
            if(makeRealPairs(PiTrig,PxRightPiClose,hDetaDphiVsPtPxRightClose)) 
                hnPiTrig->Fill(4.5);
            if(PxRightPiClose.size()!=0){ for(unsigned int i=0;i<PxRightPiClose.size();i++) hRightAssocPt->Fill(PxRightPiClose[i].Pt());}
            if(makeRealPairs(PiTrig,PxRightPiFar,hDetaDphiVsPtPxRightFar)) 
                hnPiTrig->Fill(5.5);
            if(PxRightPiFar.size()!=0){ for(unsigned int i=0;i<PxRightPiFar.size();i++) hRightAssocPt->Fill(PxRightPiFar[i].Pt());}
        }
        if(IsPxL && IsPxR){
            hnPiTrig->Fill(6.5);
        }

        if(debug)cout<<"mix event"<<endl;
        makeMixPairs(PiTrig,hDetaDphiVsPtMix);
        copyToBuffer(PiCandidate);

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

    hVertexYvsX->Fill(tpcVx,tpcVz);
    hVpdVzvsTpcVz->Fill(vpdVz,tpcVz);
    hVzDiff->Fill(vzDiff);
    hrefMult->Fill(grefMult);
    hGRefMultvsGRefMultCorr->Fill(grefMultCorr,grefMult);
    hCentrality->Fill(centrality);

    return kTRUE;
}
//-------------------------------------------
Bool_t passPion(miniDst* evt)
{
    if(debug) cout<<"IsPion loop: "<<NPiTrks<<endl;
    if(NPiTrks<=0) return kFALSE;
    for(Int_t i=0;i<NPiTrks;i++){
        Float_t pt = evt -> mPiPt[i];//trigger particle
        Float_t eta = evt -> mPiEta[i];
        Float_t phi = evt -> mPiPhi[i];
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
        Float_t eta = vec[i].Eta();
        if(IsPxL){
            if(eta<0 && eta>-0.5) PxLeftPiClose.push_back(vec[i]);
            else if(eta>0 && eta<0.5) PxLeftPiFar.push_back(vec[i]);
            if(debug)cout<<"passPionForReal IsPxL"<<endl;
        }
        if(IsPxR){
            if(eta<0 && eta>-0.5) PxRightPiFar.push_back(vec[i]);
            else if(eta>0 && eta<0.5) PxRightPiClose.push_back(vec[i]);
            if(debug)cout<<"passPionForReal IsPxR"<<endl;
        }

    }
    return kTRUE;
}
//-------------------------------------------
Int_t GetCenBinPointer(Short_t cen)
{
    //const Int_t Cen[9+1]={0,2,4,6,8,10,11,13,14,15};
    Int_t i=-999;
    if(cen>=0 && cen<=1) i=0;
    else if(cen>=2 && cen<=3) i=1; 
    else if(cen>=4 && cen<=5) i=2; 
    else if(cen>=6 && cen<=7) i=3; 
    else if(cen>=8 && cen<=9) i=4; 
    else if(cen>=10 && cen<=11) i=5; 
    else if(cen>=12 && cen<=13) i=6; 
    else if(cen==14) i=7; 
    else if(cen==15) i=8; 
    return i;
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
    hCentrality->Write();

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
