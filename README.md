# CorrelationFunction
analysis code using minitree

StMiniTreeMaker:
/star/u/liuzhen/Run14/JetLikeCorrelation/Correlation/StRoot/StMiniTreeMaker/

minitree structure:7/27/2016
Int_t    mRunId;
Int_t    mEventId;
Int_t    mGRefMult;
Float_t    mGRefMultCorr;
Float_t    mEvtWeight;
Short_t    mCentrality16;
Short_t    mCentrality9;
Float_t  mZDCRate;
Float_t  mBBCRate;
Float_t  mBField;
Float_t  mVpdVz;
Float_t  mTpcVx;
Float_t  mTpcVy;
Float_t  mTpcVz;
Float_t  mPxL; -> mPxL[mMax]
Float_t  mPxR; -> mPxR[mMax] -> every match with MTD particle have Px value

//track information
Short_t    mNTrk;
Short_t    mTrkId[mMax];
Float_t    mnSigmaPi[mMax];
Float_t    mnSigmaK[mMax];
Float_t    mnSigmaP[mMax];
Float_t    mnSigmaE[mMax];
Float_t    mgDca[mMax];
Float_t  mTrkPt[mMax];
Float_t  mTrkEta[mMax];
Float_t  mTrkPhi[mMax];


