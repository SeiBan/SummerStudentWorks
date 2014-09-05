//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  3 13:42:23 2014 by ROOT version 5.34/19
// from TTree eventsTree/new eventsTree
// found on file: muonSeedTree_1500evts.root
//////////////////////////////////////////////////////////

#ifndef MuonReconstruct_1500evts_99per_h
#define MuonReconstruct_1500evts_99per_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class MuonReconstruct_1500evts_99per {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *T_HitsAll_Muon_globalSTAseedx;
   vector<float>   *T_HitsAll_Muon_globalSTAseedy;
   vector<float>   *T_HitsAll_Muon_globalSTAseedz;
   vector<float>   *T_HitsAll_Muon_globalDirSTAseedx;
   vector<float>   *T_HitsAll_Muon_globalDirSTAseedy;
   vector<float>   *T_HitsAll_Muon_globalDirSTAseedz;
   vector<float>   *T_HitsAll_Muon_DTStation;
   vector<float>   *T_HitsAll_Muon_CSCstation;
   vector<float>   *T_Hits_Muon_globalSTAseedx;
   vector<float>   *T_Hits_Muon_globalSTAseedy;
   vector<float>   *T_Hits_Muon_globalSTAseedz;
   vector<float>   *T_Hits_Muon_globalDirSTAseedx;
   vector<float>   *T_Hits_Muon_globalDirSTAseedy;
   vector<float>   *T_Hits_Muon_globalDirSTAseedz;
   vector<float>   *T_Hits_Muon_DTStation;
   vector<float>   *T_Hits_Muon_CSCstation;
   vector<int>     *T_Seed_Muon_refFirstHit;
   vector<int>     *T_Seed_Muon_nHits;

   // List of branches
   TBranch        *b_T_HitsAll_Muon_globalSTAseedx;   //!
   TBranch        *b_T_HitsAll_Muon_globalSTAseedy;   //!
   TBranch        *b_T_HitsAll_Muon_globalSTAseedz;   //!
   TBranch        *b_T_HitsAll_Muon_globalDirSTAseedx;   //!
   TBranch        *b_T_HitsAll_Muon_globalDirSTAseedy;   //!
   TBranch        *b_T_HitsAll_Muon_globalDirSTAseedz;   //!
   TBranch        *b_T_HitsAll_Muon_DTStation;   //!
   TBranch        *b_T_HitsAll_Muon_CSCstation;   //!
   TBranch        *b_T_Hits_Muon_globalSTAseedx;   //!
   TBranch        *b_T_Hits_Muon_globalSTAseedy;   //!
   TBranch        *b_T_Hits_Muon_globalSTAseedz;   //!
   TBranch        *b_T_Hits_Muon_globalDirSTAseedx;   //!
   TBranch        *b_T_Hits_Muon_globalDirSTAseedy;   //!
   TBranch        *b_T_Hits_Muon_globalDirSTAseedz;   //!
   TBranch        *b_T_Hits_Muon_DTStation;   //!
   TBranch        *b_T_Hits_Muon_CSCstation;   //!
   TBranch        *b_T_Seed_Muon_refFirstHit;   //!
   TBranch        *b_T_Seed_Muon_nHits;   //!

   MuonReconstruct_1500evts_99per(TTree *tree=0);
   virtual ~MuonReconstruct_1500evts_99per();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MuonReconstruct_1500evts_99per_cxx
MuonReconstruct_1500evts_99per::MuonReconstruct_1500evts_99per(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muonSeedTree_1500evts.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("muonSeedTree_1500evts.root");
      }
      f->GetObject("eventsTree",tree);

   }
   Init(tree);
}

MuonReconstruct_1500evts_99per::~MuonReconstruct_1500evts_99per()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonReconstruct_1500evts_99per::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonReconstruct_1500evts_99per::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuonReconstruct_1500evts_99per::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   T_HitsAll_Muon_globalSTAseedx = 0;
   T_HitsAll_Muon_globalSTAseedy = 0;
   T_HitsAll_Muon_globalSTAseedz = 0;
   T_HitsAll_Muon_globalDirSTAseedx = 0;
   T_HitsAll_Muon_globalDirSTAseedy = 0;
   T_HitsAll_Muon_globalDirSTAseedz = 0;
   T_HitsAll_Muon_DTStation = 0;
   T_HitsAll_Muon_CSCstation = 0;
   T_Hits_Muon_globalSTAseedx = 0;
   T_Hits_Muon_globalSTAseedy = 0;
   T_Hits_Muon_globalSTAseedz = 0;
   T_Hits_Muon_globalDirSTAseedx = 0;
   T_Hits_Muon_globalDirSTAseedy = 0;
   T_Hits_Muon_globalDirSTAseedz = 0;
   T_Hits_Muon_DTStation = 0;
   T_Hits_Muon_CSCstation = 0;
   T_Seed_Muon_refFirstHit = 0;
   T_Seed_Muon_nHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("T_HitsAll_Muon_globalSTAseedx", &T_HitsAll_Muon_globalSTAseedx, &b_T_HitsAll_Muon_globalSTAseedx);
   fChain->SetBranchAddress("T_HitsAll_Muon_globalSTAseedy", &T_HitsAll_Muon_globalSTAseedy, &b_T_HitsAll_Muon_globalSTAseedy);
   fChain->SetBranchAddress("T_HitsAll_Muon_globalSTAseedz", &T_HitsAll_Muon_globalSTAseedz, &b_T_HitsAll_Muon_globalSTAseedz);
   fChain->SetBranchAddress("T_HitsAll_Muon_globalDirSTAseedx", &T_HitsAll_Muon_globalDirSTAseedx, &b_T_HitsAll_Muon_globalDirSTAseedx);
   fChain->SetBranchAddress("T_HitsAll_Muon_globalDirSTAseedy", &T_HitsAll_Muon_globalDirSTAseedy, &b_T_HitsAll_Muon_globalDirSTAseedy);
   fChain->SetBranchAddress("T_HitsAll_Muon_globalDirSTAseedz", &T_HitsAll_Muon_globalDirSTAseedz, &b_T_HitsAll_Muon_globalDirSTAseedz);
   fChain->SetBranchAddress("T_HitsAll_Muon_DTStation", &T_HitsAll_Muon_DTStation, &b_T_HitsAll_Muon_DTStation);
   fChain->SetBranchAddress("T_HitsAll_Muon_CSCstation", &T_HitsAll_Muon_CSCstation, &b_T_HitsAll_Muon_CSCstation);
   fChain->SetBranchAddress("T_Hits_Muon_globalSTAseedx", &T_Hits_Muon_globalSTAseedx, &b_T_Hits_Muon_globalSTAseedx);
   fChain->SetBranchAddress("T_Hits_Muon_globalSTAseedy", &T_Hits_Muon_globalSTAseedy, &b_T_Hits_Muon_globalSTAseedy);
   fChain->SetBranchAddress("T_Hits_Muon_globalSTAseedz", &T_Hits_Muon_globalSTAseedz, &b_T_Hits_Muon_globalSTAseedz);
   fChain->SetBranchAddress("T_Hits_Muon_globalDirSTAseedx", &T_Hits_Muon_globalDirSTAseedx, &b_T_Hits_Muon_globalDirSTAseedx);
   fChain->SetBranchAddress("T_Hits_Muon_globalDirSTAseedy", &T_Hits_Muon_globalDirSTAseedy, &b_T_Hits_Muon_globalDirSTAseedy);
   fChain->SetBranchAddress("T_Hits_Muon_globalDirSTAseedz", &T_Hits_Muon_globalDirSTAseedz, &b_T_Hits_Muon_globalDirSTAseedz);
   fChain->SetBranchAddress("T_Hits_Muon_DTStation", &T_Hits_Muon_DTStation, &b_T_Hits_Muon_DTStation);
   fChain->SetBranchAddress("T_Hits_Muon_CSCstation", &T_Hits_Muon_CSCstation, &b_T_Hits_Muon_CSCstation);
   fChain->SetBranchAddress("T_Seed_Muon_refFirstHit", &T_Seed_Muon_refFirstHit, &b_T_Seed_Muon_refFirstHit);
   fChain->SetBranchAddress("T_Seed_Muon_nHits", &T_Seed_Muon_nHits, &b_T_Seed_Muon_nHits);
   Notify();
}

Bool_t MuonReconstruct_1500evts_99per::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonReconstruct_1500evts_99per::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonReconstruct_1500evts_99per::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonReconstruct_1500evts_99per_cxx
