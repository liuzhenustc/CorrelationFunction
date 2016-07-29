#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cmath>

TFile *fin;
const Int_t nCenBin=9;
const Int_t Cen[nCenBin+1]={0,1,2,3,4,5,6,7,8,9};
//const Int_t Cen[nCenBin+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};


//const float PxLcut[nCenBin]={-2.5,-4.5,-7.5,-11.5,-16.5,-20.5,-27.5,-32.5,-37.5};//result
//const float PxRcut[nCenBin]={-2.5,-4.5,-7.5,-11.5,-17.5,-21.5,-28.5,-33.5,-40.5};

int cen9PxCut()
{
    fin = TFile::Open("forPxCut.histo.root","read");
    
    TH2D *hCenVsPxL = (TH2D*) fin -> Get("hCen9VsPxL");
    TH2D *hCenVsPxR = (TH2D*) fin -> Get("hCen9VsPxR");

    TCanvas *c2 = new TCanvas("c2","c2",0,0,1200,900);

    hCenVsPxL->Sumw2();
    hCenVsPxR->Sumw2();
    TH1D *hPxL[nCenBin];
    TH1D *hPxR[nCenBin];
    for(Int_t i=0;i<nCenBin;i++){
        Int_t binMin = hCenVsPxL->GetXaxis()->FindBin(Cen[i] + 1e-5);
        Int_t binMax = hCenVsPxL->GetXaxis()->FindBin(Cen[i+1] + 1e-5);
        hCenVsPxL->GetXaxis()->SetRange(binMin,binMax);
        hCenVsPxR->GetXaxis()->SetRange(binMin,binMax);
        hPxL[i]= (TH1D*)hCenVsPxL->ProjectionY(Form("hPxL%d",i));
        hPxL[i]->GetXaxis()->SetRangeUser(-70,0);
        hPxR[i]= (TH1D*)hCenVsPxR->ProjectionY(Form("hPxR%d",i));
        hPxR[i]->SetLineColor(2);
        hPxR[i]->GetXaxis()->SetRangeUser(-70,0);
    }
    hPxL[0]->SetTitle("centrality 70-80%");
    hPxL[1]->SetTitle("centrality 60-70%");
    hPxL[2]->SetTitle("centrality 50-60%");
    hPxL[3]->SetTitle("centrality 40-50%");
    hPxL[4]->SetTitle("centrality 30-40%");
    hPxL[5]->SetTitle("centrality 20-30%");
    hPxL[6]->SetTitle("centrality 10-20%");
    hPxL[7]->SetTitle("centrality 5-10%");
    hPxL[8]->SetTitle("centrality 0-5%");
    
    TLegend *leg;
    c2->Divide(3,3);
    for(Int_t i=0;i<nCenBin;i++){
        c2->cd(i+1);
        hPxL[i]->Draw("ep");
        hPxR[i]->Draw("epsame");
        leg = new TLegend(0.2,0.6,0.48,0.75);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->AddEntry(Form("hPxL%d",i)," -1.0<#eta<-0.5","l");
        leg->AddEntry(Form("hPxR%d",i)," 0.5<#eta<1.0","l");
        leg->Draw("same");
    }
    
    double pxLCut[nCenBin],pxRCut[nCenBin];
    Double_t pxCutPercentage = 0.1;
    for(Int_t i=0;i<nCenBin;i++){
        double intePxPlus = hPxL[i]->Integral();
        double sumPxPlus = 0;
        int jPxBin;
        for(jPxBin=1 ; sumPxPlus<pxCutPercentage*intePxPlus ; jPxBin++)
            sumPxPlus += hPxL[i]->GetBinContent(jPxBin);
        pxLCut[i] = hPxL[i]->GetBinCenter(jPxBin);
        
        sumPxPlus = 0;
        for(jPxBin=1 ; sumPxPlus<pxCutPercentage*intePxPlus ; jPxBin++)
            sumPxPlus += hPxR[i]->GetBinContent(jPxBin);
        pxRCut[i] = hPxR[i]->GetBinCenter(jPxBin);
    }
    
    cout<<"pxL cut:"<<endl;
    for(int i=0;i<nCenBin;++i)
        cout<<pxLCut[i]<<",";
    cout<<endl;
    
    cout<<"pxR cut:"<<endl;
    for(int i=0;i<nCenBin;++i)
        cout<<pxRCut[i]<<",";
    cout<<endl;



    return 0;
}
