#define MuonReconstruct_1500evts_99per_cxx
#include "MuonReconstruct_1500evts_99per.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
void MuonReconstruct_1500evts_99per::Loop(){
//      muonSeedTree_1500evts.root <- muonSeedTree6.root+muonSeedTree7.root+muonSeedTree8.root
//   In a ROOT session, you can do:
//      Root > .L MuonReconstruct_1500evts_99per.C
//      Root > MuonReconstruct_1500evts_99per t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  //  Long64_t nentries = 15;
  Long64_t nbytes = 0, nb = 0;
  
  float length(float x , float y , float z){
    return pow(x*x+y*y+z*z,0.5);
  }
  
  TFile* myFile = new TFile("MuonReconstruct_1500evts_99per.root","RECREATE");
  TTree* eventsTree = new TTree("eventsTree","tree");
  
  vector<float>* sta1_x;
  vector<float>* sta1_y;
  vector<float>* sta1_z;
  vector<float>* sta2_x;
  vector<float>* sta2_y;
  vector<float>* sta2_z;
  vector<float>* sta3_x;
  vector<float>* sta3_y;
  vector<float>* sta3_z;
  vector<float>* sta4_x;
  vector<float>* sta4_y;
  vector<float>* sta4_z;
  vector<float>* sta1dir_x;
  vector<float>* sta1dir_y;
  vector<float>* sta1dir_z;
  vector<float>* sta2dir_x;
  vector<float>* sta2dir_y;
  vector<float>* sta2dir_z;
  vector<float>* sta3dir_x;
  vector<float>* sta3dir_y;
  vector<float>* sta3dir_z;
  vector<float>* sta4dir_x;
  vector<float>* sta4dir_y;
  vector<float>* sta4dir_z;

  vector<float>* sta1temp_x;
  vector<float>* sta1temp_y;
  vector<float>* sta1temp_z;
  vector<float>* sta2temp_x;
  vector<float>* sta2temp_y;
  vector<float>* sta2temp_z;
  vector<float>* sta3temp_x;
  vector<float>* sta3temp_y;
  vector<float>* sta3temp_z;
  vector<float>* sta4temp_x;
  vector<float>* sta4temp_y;
  vector<float>* sta4temp_z;
  vector<float>* sta1dirtemp_x;
  vector<float>* sta1dirtemp_y;
  vector<float>* sta1dirtemp_z;
  vector<float>* sta2dirtemp_x;
  vector<float>* sta2dirtemp_y;
  vector<float>* sta2dirtemp_z;
  vector<float>* sta3dirtemp_x;
  vector<float>* sta3dirtemp_y;
  vector<float>* sta3dirtemp_z;
  vector<float>* sta4dirtemp_x;
  vector<float>* sta4dirtemp_y;
  vector<float>* sta4dirtemp_z;

  vector<float>* sizetemp;
  vector<float>* muontemp_x;
  vector<float>* muontemp_y;
  vector<float>* muontemp_z;
  vector<float>* muondirtemp_x;
  vector<float>* muondirtemp_y;
  vector<float>* muondirtemp_z;
  vector<int>* refFirstHit;
  vector<int>* nHits;

  eventsTree->Branch("T_HitsAll_Muon_globalSTAseedx","std::vector<float>",&T_HitsAll_Muon_globalSTAseedx);
  eventsTree->Branch("T_HitsAll_Muon_globalSTAseedy","std::vector<float>",&T_HitsAll_Muon_globalSTAseedy);
  eventsTree->Branch("T_HitsAll_Muon_globalSTAseedz","std::vector<float>",&T_HitsAll_Muon_globalSTAseedz);
  eventsTree->Branch("T_HitsAll_Muon_globalDirSTAseedx","std::vector<float>",&T_HitsAll_Muon_globalDirSTAseedx);
  eventsTree->Branch("T_HitsAll_Muon_globalDirSTAseedy","std::vector<float>",&T_HitsAll_Muon_globalDirSTAseedy);
  eventsTree->Branch("T_HitsAll_Muon_globalDirSTAseedz","std::vector<float>",&T_HitsAll_Muon_globalDirSTAseedz);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedx_ans","std::vector<float>",&T_Hits_Muon_globalSTAseedx);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedy_ans","std::vector<float>",&T_Hits_Muon_globalSTAseedy);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedz_ans","std::vector<float>",&T_Hits_Muon_globalSTAseedz);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedx_ans","std::vector<float>",&T_Hits_Muon_globalDirSTAseedx);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedy_ans","std::vector<float>",&T_Hits_Muon_globalDirSTAseedy);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedz_ans","std::vector<float>",&T_Hits_Muon_globalDirSTAseedz);
  eventsTree->Branch("T_Seed_Muon_refFirstHit_ans","std::vector<int>",&T_Seed_Muon_refFirstHit);
  eventsTree->Branch("T_Seed_Muon_nHits_ans","std::vector<int>",&T_Seed_Muon_nHits);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedx","std::vector<float>",&muontemp_x);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedy","std::vector<float>",&muontemp_y);
  eventsTree->Branch("T_Hits_Muon_globalSTAseedz","std::vector<float>",&muontemp_z);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedx","std::vector<float>",&muondirtemp_x);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedy","std::vector<float>",&muondirtemp_y);
  eventsTree->Branch("T_Hits_Muon_globalDirSTAseedz","std::vector<float>",&muondirtemp_z);
  eventsTree->Branch("T_Seed_Muon_refFirstHit","std::vector<int>",&refFirstHit);
  eventsTree->Branch("T_Seed_Muon_nHits","std::vector<int>",&nHits);

  for(Long64_t jentry = 0 ; jentry<nentries ; jentry++){
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    /*    
	  if(jentry == 38 || jentry == 41 || jentry == 43 || jentry == 44 || jentry == 80 || jentry == 82 || jentry == 133 || jentry == 140){
	  eventsTree->Fill();
	  continue;
	  }
    */

    cout<<"--------------------------------------------------------------"<<endl;        
    cout<<"                        "<<jentry<<"th event                      "<<endl; 
    cout<<"--------------------------------------------------------------"<<endl;        
    
    sta1_x = new vector<float>;
    sta1_y = new vector<float>;
    sta1_z = new vector<float>;
    sta2_x = new vector<float>;
    sta2_y = new vector<float>;
    sta2_z = new vector<float>;
    sta3_x = new vector<float>;
    sta3_y = new vector<float>;
    sta3_z = new vector<float>;
    sta4_x = new vector<float>;
    sta4_y = new vector<float>;
    sta4_z = new vector<float>;
    sta1dir_x = new vector<float>;
    sta1dir_y = new vector<float>;
    sta1dir_z = new vector<float>;
    sta2dir_x = new vector<float>;
    sta2dir_y = new vector<float>;
    sta2dir_z = new vector<float>;
    sta3dir_x = new vector<float>;
    sta3dir_y = new vector<float>;
    sta3dir_z = new vector<float>;
    sta4dir_x = new vector<float>;
    sta4dir_y = new vector<float>;
    sta4dir_z = new vector<float>;
    muontemp_x = new vector<float>;
    muontemp_y = new vector<float>;
    muontemp_z = new vector<float>;
    muondirtemp_x = new vector<float>;
    muondirtemp_y = new vector<float>;
    muondirtemp_z = new vector<float>;
    refFirstHit = new vector<int>;
    nHits = new vector<int>;
    
    /////////////////////////////////// Sorting each segment /////////////////////////////////////////
    for(unsigned int ii=0 ; ii<T_HitsAll_Muon_globalSTAseedx->size() ; ii++){ // Loop of all segments in an event
      if((T_HitsAll_Muon_DTStation->at(ii)==1) || (T_HitsAll_Muon_CSCstation->at(ii)==1)){
	sta1_x->push_back(T_HitsAll_Muon_globalSTAseedx->at(ii));
	sta1_y->push_back(T_HitsAll_Muon_globalSTAseedy->at(ii));
	sta1_z->push_back(T_HitsAll_Muon_globalSTAseedz->at(ii));
	sta1dir_x->push_back(T_HitsAll_Muon_globalDirSTAseedx->at(ii));
	sta1dir_y->push_back(T_HitsAll_Muon_globalDirSTAseedy->at(ii));
	sta1dir_z->push_back(T_HitsAll_Muon_globalDirSTAseedz->at(ii));
      } else if((T_HitsAll_Muon_DTStation->at(ii)==2) || (T_HitsAll_Muon_CSCstation->at(ii)==2)) {
	sta2_x->push_back(T_HitsAll_Muon_globalSTAseedx->at(ii));
	sta2_y->push_back(T_HitsAll_Muon_globalSTAseedy->at(ii));
	sta2_z->push_back(T_HitsAll_Muon_globalSTAseedz->at(ii));
	sta2dir_x->push_back(T_HitsAll_Muon_globalDirSTAseedx->at(ii));
	sta2dir_y->push_back(T_HitsAll_Muon_globalDirSTAseedy->at(ii));
	sta2dir_z->push_back(T_HitsAll_Muon_globalDirSTAseedz->at(ii));
      } else if((T_HitsAll_Muon_DTStation->at(ii)==3) || (T_HitsAll_Muon_CSCstation->at(ii)==3)) {
	sta3_x->push_back(T_HitsAll_Muon_globalSTAseedx->at(ii));
	sta3_y->push_back(T_HitsAll_Muon_globalSTAseedy->at(ii));
	sta3_z->push_back(T_HitsAll_Muon_globalSTAseedz->at(ii));
	sta3dir_x->push_back(T_HitsAll_Muon_globalDirSTAseedx->at(ii));
	sta3dir_y->push_back(T_HitsAll_Muon_globalDirSTAseedy->at(ii));
	sta3dir_z->push_back(T_HitsAll_Muon_globalDirSTAseedz->at(ii));
      } else if((T_HitsAll_Muon_DTStation->at(ii)==4) || (T_HitsAll_Muon_CSCstation->at(ii)==4)) {
	sta4_x->push_back(T_HitsAll_Muon_globalSTAseedx->at(ii));
	sta4_y->push_back(T_HitsAll_Muon_globalSTAseedy->at(ii));
	sta4_z->push_back(T_HitsAll_Muon_globalSTAseedz->at(ii));
	sta4dir_x->push_back(T_HitsAll_Muon_globalDirSTAseedx->at(ii));
	sta4dir_y->push_back(T_HitsAll_Muon_globalDirSTAseedy->at(ii));
	sta4dir_z->push_back(T_HitsAll_Muon_globalDirSTAseedz->at(ii));
      }
    } // Loop of all segments in an event
    //////////////////////////////////////////// Sorting each segment ///////////////////////////////////////////////       

    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    ///////////////////////////////////////////// Sarch with cone ///////////////////////////////////////////////////
    float reVec12_x[100][100];
    float reVec12_y[100][100];
    float reVec12_z[100][100];
    float inner_product12[100][100];
    float inner_productdir12[100][100];
    float reVec23_x[100][100];
    float reVec23_y[100][100];
    float reVec23_z[100][100];
    float inner_product23[100][100]; 
    float inner_productdir23[100][100];
    float reVec34_x[100][100];
    float reVec34_y[100][100];
    float reVec34_z[100][100];
    float inner_product34[100][100]; 
    float inner_productdir34[100][100];
    float reVec24_x[100][100];
    float reVec24_y[100][100];
    float reVec24_z[100][100];
    float inner_product24[100][100]; 
    
    float reVec13_x_134[100][100];
    float reVec13_y_134[100][100];
    float reVec13_z_134[100][100];
    float inner_product13_134[100][100]; 
    float inner_productdir13_134[100][100];
    float reVec34_x_134[100][100];
    float reVec34_y_134[100][100];
    float reVec34_z_134[100][100];
    float inner_product34_134[100][100]; 

    float reVec23_x_234[100][100];
    float reVec23_y_234[100][100];
    float reVec23_z_234[100][100];
    float inner_product23_234[100][100]; 
    float inner_productdir23_234[100][100];
    float reVec34_x_234[100][100];
    float reVec34_y_234[100][100];
    float reVec34_z_234[100][100];
    float inner_product34_234[100][100]; 

    float reVec14_x_14[100][100];
    float reVec14_y_14[100][100];
    float reVec14_z_14[100][100];
    float inner_product14_14[100][100]; 

    float reVec24_x_24[100][100];
    float reVec24_y_24[100][100];
    float reVec24_z_24[100][100];
    float inner_product24_24[100][100]; 

    float reVec34_x_34[100][100];
    float reVec34_y_34[100][100];
    float reVec34_z_34[100][100];
    float inner_product34_34[100][100]; 

    /////////////////////// sta1->sta2->sta3->sta4 ////////////////////////////
    for(int i1=sta1_x->size()-1 ; 0<=i1 ; i1--){ // Loop of all sta1 segments (1->2->3->4)
      unsigned int flag_1234_1 = 0;
      cout<<"aho1-1  i1="<<i1<<endl;
      sta1temp_x = new vector<float>;
      sta1temp_y = new vector<float>;
      sta1temp_z = new vector<float>;
      sta2temp_x = new vector<float>;
      sta2temp_y = new vector<float>;
      sta2temp_z = new vector<float>;
      sta3temp_x = new vector<float>;
      sta3temp_y = new vector<float>;
      sta3temp_z = new vector<float>;
      sta4temp_x = new vector<float>;
      sta4temp_y = new vector<float>;
      sta4temp_z = new vector<float>;
      sta1dirtemp_x = new vector<float>;
      sta1dirtemp_y = new vector<float>;
      sta1dirtemp_z = new vector<float>;
      sta2dirtemp_x = new vector<float>;
      sta2dirtemp_y = new vector<float>;
      sta2dirtemp_z = new vector<float>;
      sta3dirtemp_x = new vector<float>;
      sta3dirtemp_y = new vector<float>;
      sta3dirtemp_z = new vector<float>;
      sta4dirtemp_x = new vector<float>;
      sta4dirtemp_y = new vector<float>;
      sta4dirtemp_z = new vector<float>;
      sizetemp = new vector<float>;
      for(int i2=sta2_x->size()-1 ; 0<=i2 ; i2--){ // Loop of all sta2 segments (1->2->3->4)
	cout<<"aho1-2  i2="<<i2<<endl;
	reVec12_x[i1][i2] = sta2_x->at(i2)-sta1_x->at(i1);
	reVec12_y[i1][i2] = sta2_y->at(i2)-sta1_y->at(i1);
	reVec12_z[i1][i2] = sta2_z->at(i2)-sta1_z->at(i1);
	inner_product12[i1][i2] = (sta1dir_x->at(i1))*(reVec12_x[i1][i2])+(sta1dir_y->at(i1))*(reVec12_y[i1][i2])+(sta1dir_z->at(i1))*(reVec12_z[i1][i2]);
	inner_productdir12[i1][i2] = (sta1dir_x->at(i1))*(sta2dir_x->at(i2))+(sta1dir_y->at(i1))*(sta2dir_y->at(i2))+(sta1dir_z->at(i1))*(sta2dir_z->at(i2));	
	if( (inner_product12[i1][i2]>length(reVec12_x[i1][i2],reVec12_y[i1][i2],reVec12_z[i1][i2])*cos(66.5*3.14/180)) && (inner_productdir12[i1][i2]>cos(45*3.14/180)) ){
	  cout<<"cone 1->2 (1->2->3->4)"<<endl;
	  sta1temp_x->push_back(sta1_x->at(i1));
	  sta1temp_y->push_back(sta1_y->at(i1));
	  sta1temp_z->push_back(sta1_z->at(i1));
	  sta1dirtemp_x->push_back(sta1dir_x->at(i1));
	  sta1dirtemp_y->push_back(sta1dir_y->at(i1));
	  sta1dirtemp_z->push_back(sta1dir_z->at(i1));
	  sta2temp_x->push_back(sta2_x->at(i2));
	  sta2temp_y->push_back(sta2_y->at(i2));
	  sta2temp_z->push_back(sta2_z->at(i2));
	  sta2dirtemp_x->push_back(sta2dir_x->at(i2));
	  sta2dirtemp_y->push_back(sta2dir_y->at(i2));
	  sta2dirtemp_z->push_back(sta2dir_z->at(i2));
	  sizetemp->push_back(1);
	  sizetemp->push_back(1);
	  int flag_1234_3 = 0; 
	  for(int i3=sta3_x->size()-1 ; 0<=i3 ; i3--){ // Loop of all sta3 segments (1->2->3->4)
	    cout<<"aho1-3  i3="<<i3<<endl;
	    reVec23_x[i2][i3] = sta3_x->at(i3) - sta2_x->at(i2);
	    reVec23_y[i2][i3] = sta3_y->at(i3) - sta2_y->at(i2);
	    reVec23_z[i2][i3] = sta3_z->at(i3) - sta2_z->at(i2);
	    inner_product23[i2][i3] = (sta2dir_x->at(i2))*(reVec23_x[i2][i3])+(sta2dir_y->at(i2))*(reVec23_y[i2][i3])+(sta2dir_z->at(i2))*(reVec23_z[i2][i3]);
	    inner_productdir23[i2][i3] = (sta2dir_x->at(i2))*(sta3dir_x->at(i3))+(sta2dir_y->at(i2))*(sta3dir_y->at(i3))+(sta2dir_z->at(i2))*(sta3dir_z->at(i3));
	    if( (inner_product23[i2][i3]>length(reVec23_x[i2][i3],reVec23_y[i2][i3],reVec23_z[i2][i3])*cos(60*3.14/180)) && (inner_productdir23[i2][i3]>cos(50*3.14/180)) ){
	      cout<<"cone 2->3 (1->2->3->4)"<<endl;
	      sta3temp_x->push_back(sta3_x->at(i3));
	      sta3temp_y->push_back(sta3_y->at(i3));
	      sta3temp_z->push_back(sta3_z->at(i3));
	      sta3dirtemp_x->push_back(sta3dir_x->at(i3));
	      sta3dirtemp_y->push_back(sta3dir_y->at(i3));
	      sta3dirtemp_z->push_back(sta3dir_z->at(i3));
	      sizetemp->push_back(1);
	      flag_1234_3 = 1;
	      for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (1->2->3->4)
		cout<<"aho1-4  i4="<<i4<<endl;
		reVec34_x[i3][i4] = sta4_x->at(i4) - sta3_x->at(i3);
		reVec34_y[i3][i4] = sta4_y->at(i4) - sta3_y->at(i3);
		reVec34_z[i3][i4] = sta4_z->at(i4) - sta3_z->at(i3);
		inner_product34[i3][i4] = (sta3dir_x->at(i3))*(reVec34_x[i3][i4])+(sta3dir_y->at(i3))*(reVec34_y[i3][i4])+(sta3dir_z->at(i3))*(reVec34_z[i3][i4]);
		if(inner_product34[i3][i4]>length(reVec34_x[i3][i4],reVec34_y[i3][i4],reVec34_z[i3][i4])*cos(60*3.14/180)){
		  cout<<"cone 3->4 (1->2->3->4)"<<endl;
		  sta4temp_x->push_back(sta4_x->at(i4));
		  sta4temp_y->push_back(sta4_y->at(i4));
		  sta4temp_z->push_back(sta4_z->at(i4));
		  sta4dirtemp_x->push_back(sta4dir_x->at(i4));
		  sta4dirtemp_y->push_back(sta4dir_y->at(i4));
		  sta4dirtemp_z->push_back(sta4dir_z->at(i4));
		  sizetemp->push_back(1);
		  sta4_x->erase(sta4_x->begin()+i4);
		  sta4_y->erase(sta4_y->begin()+i4);
		  sta4_z->erase(sta4_z->begin()+i4);
		  sta4dir_x->erase(sta4dir_x->begin()+i4);
		  sta4dir_y->erase(sta4dir_y->begin()+i4);
		  sta4dir_z->erase(sta4dir_z->begin()+i4);
		} // if(cone sta3->sta4)
	      } // Loop of all sta4 segments (1->2->3->4)
	      sta3_x->erase(sta3_x->begin()+i3);
	      sta3_y->erase(sta3_y->begin()+i3);
	      sta3_z->erase(sta3_z->begin()+i3);
	      sta3dir_x->erase(sta3dir_x->begin()+i3);
	      sta3dir_y->erase(sta3dir_y->begin()+i3);
	      sta3dir_z->erase(sta3dir_z->begin()+i3);
	    } // if(cone sta2->sta3)
	  } // Loop of all sta3 segments (1->2->3->4)
	  if(flag_1234_3 == 0){ 
	    for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta3 segments (1->2->3->4) and 1->2->4
	      cout<<"aho1-5  i4="<<i4<<endl;
	      reVec24_x[i2][i4] = sta4_x->at(i4) - sta2_x->at(i2);
	      reVec24_y[i2][i4] = sta4_y->at(i4) - sta2_y->at(i2);
	      reVec24_z[i2][i4] = sta4_z->at(i4) - sta2_z->at(i2);
	      inner_product24[i2][i4] = (sta2dir_x->at(i2))*(reVec24_x[i2][i4])+(sta2dir_y->at(i2))*(reVec24_y[i2][i4])+(sta2dir_z->at(i2))*(reVec24_z[i2][i4]);
	      if(inner_product24[i2][i4]>length(reVec24_x[i2][i4],reVec24_y[i2][i4],reVec24_z[i2][i4])*cos(55*3.14/180)){
		cout<<"cone 2->4 (1->2->3->4)"<<endl;
		sta4temp_x->push_back(sta4_x->at(i4));
		sta4temp_y->push_back(sta4_y->at(i4));
		sta4temp_z->push_back(sta4_z->at(i4));
		sta4dirtemp_x->push_back(sta4dir_x->at(i4));
		sta4dirtemp_y->push_back(sta4dir_y->at(i4));
		sta4dirtemp_z->push_back(sta4dir_z->at(i4));
		sizetemp->push_back(1);
		sta4_x->erase(sta4_x->begin()+i4);
		sta4_y->erase(sta4_y->begin()+i4);
		sta4_z->erase(sta4_z->begin()+i4);	
		sta4dir_x->erase(sta4dir_x->begin()+i4);
		sta4dir_y->erase(sta4dir_y->begin()+i4);
		sta4dir_z->erase(sta4dir_z->begin()+i4);
	      } // if(cone sta2->sta4)
	    } // Loop of all sta3 segments (1->2->3->4) and 1->2->4
	  } // if(flag_1234_3 == 0)
	  sta2_x->erase(sta2_x->begin()+i2);
	  sta2_y->erase(sta2_y->begin()+i2);
	  sta2_z->erase(sta2_z->begin()+i2);
	  sta2dir_x->erase(sta2dir_x->begin()+i2);
	  sta2dir_y->erase(sta2dir_y->begin()+i2);
	  sta2dir_z->erase(sta2dir_z->begin()+i2);
	  flag_1234_1 = 1;
	} // if(cone sta1->sta2)
      } // Loop of all sta2 segments (1->2->3->4)
      if(flag_1234_1 == 1){
	sta1_x->erase(sta1_x->begin()+i1);
	sta1_y->erase(sta1_y->begin()+i1);
	sta1_z->erase(sta1_z->begin()+i1);
	sta1dir_x->erase(sta1dir_x->begin()+i1);
	sta1dir_y->erase(sta1dir_y->begin()+i1);
	sta1dir_z->erase(sta1dir_z->begin()+i1);
      }
      cout<<"aho1-6"<<endl;
      if(sizetemp->size() > 1){
	refFirstHit->push_back(muontemp_x->size());
	if(sta1temp_x->size()>0){
	  for(unsigned int ii1=0 ; ii1<sta1temp_x->size() ; ii1++){
	    muontemp_x->push_back(sta1temp_x->at(ii1));
	    muontemp_y->push_back(sta1temp_y->at(ii1));
	    muontemp_z->push_back(sta1temp_z->at(ii1));
	    muondirtemp_x->push_back(sta1dirtemp_x->at(ii1));
	    muondirtemp_y->push_back(sta1dirtemp_y->at(ii1));
	    muondirtemp_z->push_back(sta1dirtemp_z->at(ii1));
	  }
	}
	if(sta2temp_x->size()>0){
	  for(unsigned int ii2=0 ; ii2<sta2temp_x->size() ; ii2++){
	    muontemp_x->push_back(sta2temp_x->at(ii2));
	    muontemp_y->push_back(sta2temp_y->at(ii2));
	    muontemp_z->push_back(sta2temp_z->at(ii2));
	    muondirtemp_x->push_back(sta2dirtemp_x->at(ii2));
	    muondirtemp_y->push_back(sta2dirtemp_y->at(ii2));
	    muondirtemp_z->push_back(sta2dirtemp_z->at(ii2));
	  }
	}
	if(sta3temp_x->size()>0){
	  for(unsigned int ii3=0 ; ii3<sta3temp_x->size() ; ii3++){
	    muontemp_x->push_back(sta3temp_x->at(ii3));
	    muontemp_y->push_back(sta3temp_y->at(ii3));
	    muontemp_z->push_back(sta3temp_z->at(ii3));
	    muondirtemp_x->push_back(sta3dirtemp_x->at(ii3));
	    muondirtemp_y->push_back(sta3dirtemp_y->at(ii3));
	    muondirtemp_z->push_back(sta3dirtemp_z->at(ii3));
	  }
	}
	if(sta4temp_x->size()>0){
	  for(unsigned int ii4=0 ; ii4<sta4temp_x->size() ; ii4++){
	    muontemp_x->push_back(sta4temp_x->at(ii4));
	    muontemp_y->push_back(sta4temp_y->at(ii4));
	    muontemp_z->push_back(sta4temp_z->at(ii4));
	    muondirtemp_x->push_back(sta4dirtemp_x->at(ii4));
	    muondirtemp_y->push_back(sta4dirtemp_y->at(ii4));
	    muondirtemp_z->push_back(sta4dirtemp_z->at(ii4));
	  }
	}
	nHits->push_back(sizetemp->size());
      }
      delete sta1temp_x;
      delete sta1temp_y;
      delete sta1temp_z;
      delete sta2temp_x;
      delete sta2temp_y;
      delete sta2temp_z;
      delete sta3temp_x;
      delete sta3temp_y;
      delete sta3temp_z;
      delete sta4temp_x;
      delete sta4temp_y;
      delete sta4temp_z;
      delete sta1dirtemp_x;
      delete sta1dirtemp_y;
      delete sta1dirtemp_z;
      delete sta2dirtemp_x;
      delete sta2dirtemp_y;
      delete sta2dirtemp_z;
      delete sta3dirtemp_x;
      delete sta3dirtemp_y;
      delete sta3dirtemp_z;
      delete sta4dirtemp_x;
      delete sta4dirtemp_y;
      delete sta4dirtemp_z;
      delete sizetemp;
    } // Loop of all sta1 segments (1->2->3->4)
    /////////////////////// sta1->sta2->sta3->sta4 ///////////////////////////////
    
    cout<<"-------------------------After 1->2->3->4-------------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    ///////////////////////////// sta1->sta3->sta4 ////////////////////////////////////
    unsigned int check_134_1 = 1;
    unsigned int check_134_2 = 1;
    unsigned int check_134_3 = 1;
    for(unsigned int l1=0 ; l1<sta1_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta1_x->at(l1) == muontemp_x->at(l2)){
	  check_134_1 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta3_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta3_x->at(l1) == muontemp_x->at(l2)){
	  check_134_2 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta4_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta4_x->at(l1) == muontemp_x->at(l2)){
	  check_134_3 = 2;
	}
      }
    }
    cout<<"aho2-1"<<endl;
    if(check_134_1 == 1 && check_134_2 == 1 && check_134_3 == 1){
      for(int i1=sta1_x->size()-1 ; 0<=i1 ; i1--){ // Loop of all sta1 segments (1->3->4)
	cout<<"aho2-2 i1="<<i1<<endl;
	sta1temp_x = new vector<float>;
	sta1temp_y = new vector<float>;
	sta1temp_z = new vector<float>;
	sta3temp_x = new vector<float>;
	sta3temp_y = new vector<float>;
	sta3temp_z = new vector<float>;
	sta4temp_x = new vector<float>;
	sta4temp_y = new vector<float>;
	sta4temp_z = new vector<float>;
	sta1dirtemp_x = new vector<float>;
	sta1dirtemp_y = new vector<float>;
	sta1dirtemp_z = new vector<float>;
	sta3dirtemp_x = new vector<float>;
	sta3dirtemp_y = new vector<float>;
	sta3dirtemp_z = new vector<float>;
	sta4dirtemp_x = new vector<float>;
	sta4dirtemp_y = new vector<float>;
	sta4dirtemp_z = new vector<float>;
	sizetemp = new vector<float>;
	unsigned int flag_134_1 = 0;
	for(int i3=sta3_x->size()-1 ; 0<=i3 ; i3--){ // Loop of all sta3 segments (1->3->4)
	  cout<<"aho2-3 i3="<<i3<<endl;
	  reVec13_x_134[i1][i3] = sta3_x->at(i3)-sta1_x->at(i1);
	  reVec13_y_134[i1][i3] = sta3_y->at(i3)-sta1_y->at(i1);
	  reVec13_z_134[i1][i3] = sta3_z->at(i3)-sta1_z->at(i1);
	  inner_product13_134[i1][i3] = (sta1dir_x->at(i1))*(reVec13_x_134[i1][i3])+(sta1dir_y->at(i1))*(reVec13_y_134[i1][i3])+(sta1dir_z->at(i1))*(reVec13_z_134[i1][i3]);
	  inner_productdir13_134[i1][i3] = (sta1dir_x->at(i1))*(sta3dir_x->at(i3))+(sta1dir_y->at(i1))*(sta3dir_y->at(i3))+(sta1dir_z->at(i1))*(sta3dir_z->at(i3));	
	  if( (inner_product13_134[i1][i3]>length(reVec13_x_134[i1][i3],reVec13_y_134[i1][i3],reVec13_z_134[i1][i3])*cos(54*3.14/180)) && (inner_productdir13_134[i1][i3])>cos(50*3.14/180) ){
	    cout<<"cone 1->3 (1->3->4)"<<endl;
	    sta1temp_x->push_back(sta1_x->at(i1));
	    sta1temp_y->push_back(sta1_y->at(i1));
	    sta1temp_z->push_back(sta1_z->at(i1));
	    sta1dirtemp_x->push_back(sta1dir_x->at(i1));
	    sta1dirtemp_y->push_back(sta1dir_y->at(i1));
	    sta1dirtemp_z->push_back(sta1dir_z->at(i1));
	    sta3temp_x->push_back(sta3_x->at(i3));
	    sta3temp_y->push_back(sta3_y->at(i3));
	    sta3temp_z->push_back(sta3_z->at(i3));
	    sta3dirtemp_x->push_back(sta3dir_x->at(i3));
	    sta3dirtemp_y->push_back(sta3dir_y->at(i3));
	    sta3dirtemp_z->push_back(sta3dir_z->at(i3));
	    sizetemp->push_back(1);
	    sizetemp->push_back(1);
	    for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (1->3->4)
	      cout<<"aho2-4 i4="<<i4<<endl;
	      reVec34_x_134[i3][i4] = sta4_x->at(i4) - sta3_x->at(i3);
	      reVec34_y_134[i3][i4] = sta4_y->at(i4) - sta3_y->at(i3);
	      reVec34_z_134[i3][i4] = sta4_z->at(i4) - sta3_z->at(i3);
	      inner_product34_134[i3][i4] = (sta3dir_x->at(i3))*(reVec34_x_134[i3][i4])+(sta3dir_y->at(i3))*(reVec34_y_134[i3][i4])+(sta3dir_z->at(i3))*(reVec34_z_134[i3][i4]);
	      if( (inner_product34_134[i3][i4]>length(reVec34_x_134[i3][i4],reVec34_y_134[i3][i4],reVec34_z_134[i3][i4])*cos(60*3.14/180)) ){
		cout<<"cone 3->4 (1->3->4)"<<endl;
		sta4temp_x->push_back(sta4_x->at(i4));
		sta4temp_y->push_back(sta4_y->at(i4));
		sta4temp_z->push_back(sta4_z->at(i4));
		sta4dirtemp_x->push_back(sta4dir_x->at(i4));
		sta4dirtemp_y->push_back(sta4dir_y->at(i4));
		sta4dirtemp_z->push_back(sta4dir_z->at(i4));
		sizetemp->push_back(1);
		sta4_x->erase(sta4_x->begin()+i4);
		sta4_y->erase(sta4_y->begin()+i4);
		sta4_z->erase(sta4_z->begin()+i4);
		sta4dir_x->erase(sta4dir_x->begin()+i4);
		sta4dir_y->erase(sta4dir_y->begin()+i4);
		sta4dir_z->erase(sta4dir_z->begin()+i4);
	      }
	    } // Loop of all sta4 segments (1->3->4)
	    sta3_x->erase(sta3_x->begin()+i3);
	    sta3_y->erase(sta3_y->begin()+i3);
	    sta3_z->erase(sta3_z->begin()+i3);
	    sta3dir_x->erase(sta3dir_x->begin()+i3);
	    sta3dir_y->erase(sta3dir_y->begin()+i3);
	    sta3dir_z->erase(sta3dir_z->begin()+i3);
	    flag_134_1 = 1;
	  } // if(cone search 1->3)
	} // Loop of all sta3 segments (1->3->4)
	if(flag_134_1 == 1){
	  sta1_x->erase(sta1_x->begin()+i1);
	  sta1_y->erase(sta1_y->begin()+i1);
	  sta1_z->erase(sta1_z->begin()+i1);
	  sta1dir_x->erase(sta1dir_x->begin()+i1);
	  sta1dir_y->erase(sta1dir_y->begin()+i1);
	  sta1dir_z->erase(sta1dir_z->begin()+i1);
	}
	cout<<"aho2-5"<<endl;
	if(sizetemp->size() > 1){
	  refFirstHit->push_back(muontemp_x->size());
	  if(sta1temp_x->size()>0){
	    for(unsigned int ii1=0 ; ii1<sta1temp_x->size() ; ii1++){
	      muontemp_x->push_back(sta1temp_x->at(ii1));
	      muontemp_y->push_back(sta1temp_y->at(ii1));
	      muontemp_z->push_back(sta1temp_z->at(ii1));
	      muondirtemp_x->push_back(sta1dirtemp_x->at(ii1));
	      muondirtemp_y->push_back(sta1dirtemp_y->at(ii1));
	      muondirtemp_z->push_back(sta1dirtemp_z->at(ii1));
	    }
	  }
	  if(sta3temp_x->size()>0){
	    for(unsigned int ii3=0 ; ii3<sta3temp_x->size() ; ii3++){
	      muontemp_x->push_back(sta3temp_x->at(ii3));
	      muontemp_y->push_back(sta3temp_y->at(ii3));
	      muontemp_z->push_back(sta3temp_z->at(ii3));
	      muondirtemp_x->push_back(sta3dirtemp_x->at(ii3));
	      muondirtemp_y->push_back(sta3dirtemp_y->at(ii3));
	      muondirtemp_z->push_back(sta3dirtemp_z->at(ii3));
	    }
	  }
	  if(sta4temp_x->size()>0){
	    for(unsigned int ii4=0 ; ii4<sta4temp_x->size() ; ii4++){
	      muontemp_x->push_back(sta4temp_x->at(ii4));
	      muontemp_y->push_back(sta4temp_y->at(ii4));
	      muontemp_z->push_back(sta4temp_z->at(ii4));
	      muondirtemp_x->push_back(sta4dirtemp_x->at(ii4));
	      muondirtemp_y->push_back(sta4dirtemp_y->at(ii4));
	      muondirtemp_z->push_back(sta4dirtemp_z->at(ii4));
	    }
	  }
	  nHits->push_back(sizetemp->size());
	}
	delete sta1temp_x;
	delete sta1temp_y;
	delete sta1temp_z;
	delete sta3temp_x;
	delete sta3temp_y;
	delete sta3temp_z;
	delete sta4temp_x;
	delete sta4temp_y;
	delete sta4temp_z;
	delete sta1dirtemp_x;
	delete sta1dirtemp_y;
	delete sta1dirtemp_z;
	delete sta3dirtemp_x;
	delete sta3dirtemp_y;
	delete sta3dirtemp_z;
	delete sta4dirtemp_x;
	delete sta4dirtemp_y;
	delete sta4dirtemp_z;
	delete sizetemp;
      } // Loop of all sta1 segments (1->3->4)
    }
    cout<<"aho2-6"<<endl;
    ///////////////////////sta1->sta3->sta4////////////////////////////////////

    cout<<"--------------------------After 1->3->4------------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    /////////////////////////// sta2->sta3->sta4 /////////////////////////////
    unsigned int check_234_1 = 1;
    unsigned int check_234_2 = 1;
    unsigned int check_234_3 = 1;
    for(unsigned int l1=0 ; l1<sta2_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta2_x->at(l1) == muontemp_x->at(l2)){
	  check_234_1 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta3_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta3_x->at(l1) == muontemp_x->at(l2)){
	  check_234_2 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta4_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta4_x->at(l1) == muontemp_x->at(l2)){
	  check_234_3 = 2;
	}
      }
    }
    cout<<"aho3-1"<<endl;
    if(check_234_1 == 1 && check_234_2 == 1 && check_234_3 == 1){
      for(int i2=sta2_x->size()-1 ; 0<=i2 ; i2--){ // Loop of all sta2 segments (2->3->4)
	cout<<"aho3-2 i2="<<i2<<endl;
	sta2temp_x = new vector<float>;
	sta2temp_y = new vector<float>;
	sta2temp_z = new vector<float>;
	sta3temp_x = new vector<float>;
	sta3temp_y = new vector<float>;
	sta3temp_z = new vector<float>;
	sta4temp_x = new vector<float>;
	sta4temp_y = new vector<float>;
	sta4temp_z = new vector<float>;
	sta2dirtemp_x = new vector<float>;
	sta2dirtemp_y = new vector<float>;
	sta2dirtemp_z = new vector<float>;
	sta3dirtemp_x = new vector<float>;
	sta3dirtemp_y = new vector<float>;
	sta3dirtemp_z = new vector<float>;
	sta4dirtemp_x = new vector<float>;
	sta4dirtemp_y = new vector<float>;
	sta4dirtemp_z = new vector<float>;
	sizetemp = new vector<float>;
	unsigned int flag_234_1 = 0;
	for(int i3=sta3_x->size()-1 ; 0<=i3 ; i3--){ // Loop of all sta3 segments (2->3->4)
	  cout<<"aho3-3 i3="<<i3<<endl;
	  reVec23_x_234[i2][i3] = sta3_x->at(i3) - sta2_x->at(i2);
	  reVec23_y_234[i2][i3] = sta3_y->at(i3) - sta2_y->at(i2);
	  reVec23_z_234[i2][i3] = sta3_z->at(i3) - sta2_z->at(i2);
	  inner_product23_234[i2][i3] = (sta2dir_x->at(i2))*(reVec23_x_234[i2][i3])+(sta2dir_y->at(i2))*(reVec23_y_234[i2][i3])+(sta2dir_z->at(i2))*(reVec23_z_234[i2][i3]);
	  inner_productdir23_234[i2][i3] = (sta2dir_x->at(i2))*(sta3dir_x->at(i3))+(sta2dir_y->at(i2))*(sta3dir_y->at(i3))+(sta2dir_z->at(i2))*(sta3dir_z->at(i3));
	  if( (inner_product23_234[i2][i3]>length(reVec23_x_234[i2][i3],reVec23_y_234[i2][i3],reVec23_z_234[i2][i3])*cos(60*3.14/180)) && (inner_productdir23_234[i2][i3])>cos(50*3.14/180) ){
	    cout<<"cone 2->3 (2->3->4)"<<endl;
	    sta2temp_x->push_back(sta2_x->at(i2));
	    sta2temp_y->push_back(sta2_y->at(i2));
	    sta2temp_z->push_back(sta2_z->at(i2));
	    sta2dirtemp_x->push_back(sta2dir_x->at(i2));
	    sta2dirtemp_y->push_back(sta2dir_y->at(i2));
	    sta2dirtemp_z->push_back(sta2dir_z->at(i2));
	    sta3temp_x->push_back(sta3_x->at(i3));
	    sta3temp_y->push_back(sta3_y->at(i3));
	    sta3temp_z->push_back(sta3_z->at(i3));
	    sta3dirtemp_x->push_back(sta3dir_x->at(i3));
	    sta3dirtemp_y->push_back(sta3dir_y->at(i3));
	    sta3dirtemp_z->push_back(sta3dir_z->at(i3));
	    sizetemp->push_back(1);
	    sizetemp->push_back(1);
	    for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (2->3->4)
	      cout<<"aho2-4 i4="<<i4<<endl;
	      reVec34_x_234[i3][i4] = sta4_x->at(i4) - sta3_x->at(i3);
	      reVec34_y_234[i3][i4] = sta4_y->at(i4) - sta3_y->at(i3);
	      reVec34_z_234[i3][i4] = sta4_z->at(i4) - sta3_z->at(i3);
	      inner_product34_234[i3][i4] = (sta2dir_x->at(i2))*(reVec34_x_234[i3][i4])+(sta2dir_y->at(i2))*(reVec34_y_234[i3][i4])+(sta2dir_z->at(i2))*(reVec34_z_234[i3][i4]);
	      if( (inner_product34_234[i3][i4]>length(reVec34_x_234[i3][i4],reVec34_y_234[i3][i4],reVec34_z_234[i3][i4])*cos(60*3.14/180)) ){
		cout<<"cone 3->4 (2->3->4)"<<endl;
		sta4temp_x->push_back(sta4_x->at(i4));
		sta4temp_y->push_back(sta4_y->at(i4));
		sta4temp_z->push_back(sta4_z->at(i4));
		sta4dirtemp_x->push_back(sta4dir_x->at(i4));
		sta4dirtemp_y->push_back(sta4dir_y->at(i4));
		sta4dirtemp_z->push_back(sta4dir_z->at(i4));
		sizetemp->push_back(1);
		sta4_x->erase(sta4_x->begin()+i4);
		sta4_y->erase(sta4_y->begin()+i4);
		sta4_z->erase(sta4_z->begin()+i4);
		sta4dir_x->erase(sta4dir_x->begin()+i4);
		sta4dir_y->erase(sta4dir_y->begin()+i4);
		sta4dir_z->erase(sta4dir_z->begin()+i4);
	      }
	    } // Loop of all sta4 segments (2->3->4)
	    sta3_x->erase(sta3_x->begin()+i3);
	    sta3_y->erase(sta3_y->begin()+i3);
	    sta3_z->erase(sta3_z->begin()+i3);
	    sta3dir_x->erase(sta3dir_x->begin()+i3);
	    sta3dir_y->erase(sta3dir_y->begin()+i3);
	    sta3dir_z->erase(sta3dir_z->begin()+i3);
	    flag_234_1 = 1;
	  }
	} // Loop of all sta3 segments (2->3->4)
	if(flag_234_1 == 1){
	  sta2_x->erase(sta2_x->begin()+i2);
	  sta2_y->erase(sta2_y->begin()+i2);
	  sta2_z->erase(sta2_z->begin()+i2);
	  sta2dir_x->erase(sta2dir_x->begin()+i2);
	  sta2dir_y->erase(sta2dir_y->begin()+i2);
	  sta2dir_z->erase(sta2dir_z->begin()+i2);
	}
	cout<<"aho3-5"<<endl;
	if(sizetemp->size()>1){
	  refFirstHit->push_back(muontemp_x->size());
	  if(sta2temp_x->size()>0){
	    for(unsigned int ii2=0 ; ii2<sta2temp_x->size() ; ii2++){
	      muontemp_x->push_back(sta2temp_x->at(ii2));
	      muontemp_y->push_back(sta2temp_y->at(ii2));
	      muontemp_z->push_back(sta2temp_z->at(ii2));
	      muondirtemp_x->push_back(sta2dirtemp_x->at(ii2));
	      muondirtemp_y->push_back(sta2dirtemp_y->at(ii2));
	      muondirtemp_z->push_back(sta2dirtemp_z->at(ii2));
	    }
	  }
	  if(sta3temp_x->size()>0){
	    for(unsigned int ii3=0 ; ii3<sta3temp_x->size() ; ii3++){
	      muontemp_x->push_back(sta3temp_x->at(ii3));
	      muontemp_y->push_back(sta3temp_y->at(ii3));
	      muontemp_z->push_back(sta3temp_z->at(ii3));
	      muondirtemp_x->push_back(sta3dirtemp_x->at(ii3));
	      muondirtemp_y->push_back(sta3dirtemp_y->at(ii3));
	      muondirtemp_z->push_back(sta3dirtemp_z->at(ii3));
	    }
	  }
	  if(sta4temp_x->size()>0){
	    for(unsigned int ii4=0 ; ii4<sta4temp_x->size() ; ii4++){
	      muontemp_x->push_back(sta4temp_x->at(ii4));
	      muontemp_y->push_back(sta4temp_y->at(ii4));
	      muontemp_z->push_back(sta4temp_z->at(ii4));
	      muondirtemp_x->push_back(sta4dirtemp_x->at(ii4));
	      muondirtemp_y->push_back(sta4dirtemp_y->at(ii4));
	      muondirtemp_z->push_back(sta4dirtemp_z->at(ii4));
	    }
	  }
	  nHits->push_back(sizetemp->size());
	}
	delete sta2temp_x;
	delete sta2temp_y;
	delete sta2temp_z;
	delete sta3temp_x;
	delete sta3temp_y;
	delete sta3temp_z;
	delete sta4temp_x;
	delete sta4temp_y;
	delete sta4temp_z;
	delete sta2dirtemp_x;
	delete sta2dirtemp_y;
	delete sta2dirtemp_z;
	delete sta3dirtemp_x;
	delete sta3dirtemp_y;
	delete sta3dirtemp_z;
	delete sta4dirtemp_x;
	delete sta4dirtemp_y;
	delete sta4dirtemp_z;
	delete sizetemp;
      } // Loop of all sta2 segments (2->3->4)
    }
    cout<<"aho3-6"<<endl;
    /////////////////////////// sta2->sta3->sta4 /////////////////////////////    
    
    cout<<"----------------------------After 2->3->4----------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    /////////////////////////// sta1->sta4 /////////////////////////////
    unsigned int check_14_1 = 1;
    unsigned int check_14_2 = 1;
    for(unsigned int l1=0 ; l1<sta1_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta1_x->at(l1) == muontemp_x->at(l2)){
	  check_14_1 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta4_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta4_x->at(l1) == muontemp_x->at(l2)){
	  check_14_2 = 2;
	}
      }
    }
    if(check_14_1 == 1 && check_14_2 == 1){
      for(int i1=sta1_x->size()-1 ; 0<=i1 ; i1--){ // Loop of all sta1 segments (1->4)
	sta1temp_x = new vector<float>;
	sta1temp_y = new vector<float>;
	sta1temp_z = new vector<float>;
	sta1dirtemp_x = new vector<float>;
	sta1dirtemp_y = new vector<float>;
	sta1dirtemp_z = new vector<float>;
	sta4temp_x = new vector<float>;
	sta4temp_y = new vector<float>;
	sta4temp_z = new vector<float>;
	sta4dirtemp_x = new vector<float>;
	sta4dirtemp_y = new vector<float>;
	sta4dirtemp_z = new vector<float>;
	sizetemp = new vector<float>;
	int flag_14_1 = 0;
	for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (1->4)
	  reVec14_x_14[i1][i4] = sta4_x->at(i4) - sta1_x->at(i1);
	  reVec14_y_14[i1][i4] = sta4_y->at(i4) - sta1_y->at(i1);
	  reVec14_z_14[i1][i4] = sta4_z->at(i4) - sta1_z->at(i1);
	  inner_product14_14[i1][i4] = (sta1dir_x->at(i1))*(reVec14_x_14[i1][i4])+(sta1dir_y->at(i1))*(reVec14_y_14[i1][i4])+(sta1dir_z->at(i1))*(reVec14_z_14[i1][i4]);
	  if( (inner_product14_14[i1][i4]>length(reVec14_x_14[i1][i4],reVec14_y_14[i1][i4],reVec14_z_14[i1][i4])*cos(50*3.14/180)) ){
	    cout<<"cone 1->4 (1->4)"<<endl;
	    flag_14_1 = 1;
	    sta1temp_x->push_back(sta1_x->at(i1));
	    sta1temp_y->push_back(sta1_y->at(i1));
	    sta1temp_z->push_back(sta1_z->at(i1));
	    sta1dirtemp_x->push_back(sta1dir_x->at(i1));
	    sta1dirtemp_y->push_back(sta1dir_y->at(i1));
	    sta1dirtemp_z->push_back(sta1dir_z->at(i1));
	    sta4temp_x->push_back(sta4_x->at(i4));
	    sta4temp_y->push_back(sta4_y->at(i4));
	    sta4temp_z->push_back(sta4_z->at(i4));
	    sta4dirtemp_x->push_back(sta4dir_x->at(i4));
	    sta4dirtemp_y->push_back(sta4dir_y->at(i4));
	    sta4dirtemp_z->push_back(sta4dir_z->at(i4));
	    sizetemp->push_back(1);
	    sizetemp->push_back(1);
	    sta4_x->erase(sta4_x->begin()+i4);
	    sta4_y->erase(sta4_y->begin()+i4);
	    sta4_z->erase(sta4_z->begin()+i4);
	    sta4dir_x->erase(sta4dir_x->begin()+i4);
	    sta4dir_y->erase(sta4dir_y->begin()+i4);
	    sta4dir_z->erase(sta4dir_z->begin()+i4);
	  }
	} // Loop of all sta4 segments (1->4)
	if(flag_14_1 == 1){
	  sta1_x->erase(sta1_x->begin()+i1);
	  sta1_y->erase(sta1_y->begin()+i1);
	  sta1_z->erase(sta1_z->begin()+i1);
	  sta1dir_x->erase(sta1dir_x->begin()+i1);
	  sta1dir_y->erase(sta1dir_y->begin()+i1);
	  sta1dir_z->erase(sta1dir_z->begin()+i1);	
	}
	if(sizetemp->size()>1){
	  refFirstHit->push_back(muontemp_x->size());
	  if(sta1temp_x->size()>0){
	    for(unsigned int ii11=0 ; ii11<sta1temp_x->size() ; ii11++){
	      muontemp_x->push_back(sta1temp_x->at(ii11));
	      muontemp_y->push_back(sta1temp_y->at(ii11));
	      muontemp_z->push_back(sta1temp_z->at(ii11));
	      muondirtemp_x->push_back(sta1dirtemp_x->at(ii11));
	      muondirtemp_y->push_back(sta1dirtemp_y->at(ii11));
	      muondirtemp_z->push_back(sta1dirtemp_z->at(ii11));
	    }
	  }
	  if(sta4temp_x->size()>0){
	    for(unsigned int ii44=0 ; ii44<sta4temp_x->size() ; ii44++){
	      muontemp_x->push_back(sta4temp_x->at(ii44));
	      muontemp_y->push_back(sta4temp_y->at(ii44));
	      muontemp_z->push_back(sta4temp_z->at(ii44));
	      muondirtemp_x->push_back(sta4dirtemp_x->at(ii44));
	      muondirtemp_y->push_back(sta4dirtemp_y->at(ii44));
	      muondirtemp_z->push_back(sta4dirtemp_z->at(ii44));
	    }
	  }
	  nHits->push_back(sizetemp->size());
	}
	delete sta1temp_x;
	delete sta1temp_y;
	delete sta1temp_z;
	delete sta4temp_x;
	delete sta4temp_y;
	delete sta4temp_z;
	delete sta1dirtemp_x;
	delete sta1dirtemp_y;
	delete sta1dirtemp_z;
	delete sta4dirtemp_x;
	delete sta4dirtemp_y;
	delete sta4dirtemp_z;
	delete sizetemp;
      } // Loop of all sta1 segments (1->4)
    }
    ////////////////////// sta1->sta4 ////////////////////////////////////// 
    cout<<"---------------------------After 1->4-----------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;

    /////////////////////////// sta2->sta4 /////////////////////////////
    unsigned int check_24_1 = 1;
    unsigned int check_24_2 = 1;
    for(unsigned int l1=0 ; l1<sta2_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta2_x->at(l1) == muontemp_x->at(l2)){
	  check_24_1 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta4_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta4_x->at(l1) == muontemp_x->at(l2)){
	  check_24_2 = 2;
	}
      }
    }
    if(check_24_1 == 1 && check_24_2 == 1){
      for(int i2=sta2_x->size()-1 ; 0<=i2 ; i2--){ // Loop of all sta2 segments (2->4)
	sta2temp_x = new vector<float>;
	sta2temp_y = new vector<float>;
	sta2temp_z = new vector<float>;
	sta4temp_x = new vector<float>;
	sta4temp_y = new vector<float>;
	sta4temp_z = new vector<float>;
	sta2dirtemp_x = new vector<float>;
	sta2dirtemp_y = new vector<float>;
	sta2dirtemp_z = new vector<float>;
	sta4dirtemp_x = new vector<float>;
	sta4dirtemp_y = new vector<float>;
	sta4dirtemp_z = new vector<float>;
	sizetemp = new vector<float>;
	int flag_24_1 = 0;
	for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (2->4)
	  reVec24_x_24[i2][i4] = sta4_x->at(i4) - sta2_x->at(i2);
	  reVec24_y_24[i2][i4] = sta4_y->at(i4) - sta2_y->at(i2);
	  reVec24_z_24[i2][i4] = sta4_z->at(i4) - sta2_z->at(i2);
	  inner_product24_24[i2][i4] = (sta2dir_x->at(i2))*(reVec24_x_24[i2][i4])+(sta2dir_y->at(i2))*(reVec24_y_24[i2][i4])+(sta2dir_z->at(i2))*(reVec24_z_24[i2][i4]);
	  if( (inner_product24_24[i2][i4]>length(reVec24_x_24[i2][i4],reVec24_y_24[i2][i4],reVec24_z_24[i2][i4])*cos(55*3.14/180)) ){
	    cout<<"cone 2->4 (2->4)"<<endl;
	    flag_24_1 = 1;
	    sta2temp_x->push_back(sta2_x->at(i2));
	    sta2temp_y->push_back(sta2_y->at(i2));
	    sta2temp_z->push_back(sta2_z->at(i2));
	    sta4temp_x->push_back(sta4_x->at(i4));
	    sta4temp_y->push_back(sta4_y->at(i4));
	    sta4temp_z->push_back(sta4_z->at(i4));
	    sta2dirtemp_x->push_back(sta2dir_x->at(i2));
	    sta2dirtemp_y->push_back(sta2dir_y->at(i2));
	    sta2dirtemp_z->push_back(sta2dir_z->at(i2));
	    sta4dirtemp_x->push_back(sta4dir_x->at(i4));
	    sta4dirtemp_y->push_back(sta4dir_y->at(i4));
	    sta4dirtemp_z->push_back(sta4dir_z->at(i4));
	    sizetemp->push_back(1);
	    sizetemp->push_back(1);
	    sta4_x->erase(sta4_x->begin()+i4);
	    sta4_y->erase(sta4_y->begin()+i4);
	    sta4_z->erase(sta4_z->begin()+i4);
	    sta4dir_x->erase(sta4dir_x->begin()+i4);
	    sta4dir_y->erase(sta4dir_y->begin()+i4);
	    sta4dir_z->erase(sta4dir_z->begin()+i4);
	  }
	} // Loop of all sta4 segments (2->4)
	if(flag_24_1 == 1){
	  sta2_x->erase(sta2_x->begin()+i2);
	  sta2_y->erase(sta2_y->begin()+i2);
	  sta2_z->erase(sta2_z->begin()+i2);
	  sta2dir_x->erase(sta2dir_x->begin()+i2);
	  sta2dir_y->erase(sta2dir_y->begin()+i2);
	  sta2dir_z->erase(sta2dir_z->begin()+i2);	
	}
	if(sizetemp->size()>1){
	  refFirstHit->push_back(muontemp_x->size());
	  if(sta2temp_x->size()>0){
	    for(unsigned int ii22=0 ; ii22<sta2temp_x->size() ; ii22++){
	      muontemp_x->push_back(sta2temp_x->at(ii22));
	      muontemp_y->push_back(sta2temp_y->at(ii22));
	      muontemp_z->push_back(sta2temp_z->at(ii22));
	      muondirtemp_x->push_back(sta2dirtemp_x->at(ii22));
	      muondirtemp_y->push_back(sta2dirtemp_y->at(ii22));
	      muondirtemp_z->push_back(sta2dirtemp_z->at(ii22));
	    }
	  }
	  if(sta4temp_x->size()>0){
	    for(unsigned int ii44=0 ; ii44<sta4temp_x->size() ; ii44++){
	      muontemp_x->push_back(sta4temp_x->at(ii44));
	      muontemp_y->push_back(sta4temp_y->at(ii44));
	      muontemp_z->push_back(sta4temp_z->at(ii44));
	      muondirtemp_x->push_back(sta4dirtemp_x->at(ii44));
	      muondirtemp_y->push_back(sta4dirtemp_y->at(ii44));
	      muondirtemp_z->push_back(sta4dirtemp_z->at(ii44));
	    }
	  }
	  nHits->push_back(sizetemp->size());
	}
	delete sta2temp_x;
	delete sta2temp_y;
	delete sta2temp_z;
	delete sta4temp_x;
	delete sta4temp_y;
	delete sta4temp_z;
	delete sta2dirtemp_x;
	delete sta2dirtemp_y;
	delete sta2dirtemp_z;
	delete sta4dirtemp_x;
	delete sta4dirtemp_y;
	delete sta4dirtemp_z;
	delete sizetemp;
      } // Loop of all sta2 segments (2->4)
    }
    ////////////////////// sta2->sta4 ////////////////////////////////////// 
    cout<<"---------------------------After 2->4-----------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    /////////////////////////// sta3->sta4 /////////////////////////////
    unsigned int check_34_1 = 1;
    unsigned int check_34_2 = 1;
    for(unsigned int l1=0 ; l1<sta3_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta3_x->at(l1) == muontemp_x->at(l2)){
	  check_34_1 = 2;
	}
      }
    }
    for(unsigned int l1=0 ; l1<sta4_x->size() ; l1++){
      for(unsigned int l2=0 ; l2<muontemp_x->size() ; l2++){
	if(sta4_x->at(l1) == muontemp_x->at(l2)){
	  check_34_2 = 2;
	}
      }
    }
    if(check_34_1 == 1 && check_34_2 == 1){
      for(int i3=sta3_x->size()-1 ; 0<=i3 ; i3--){ // Loop of all sta3 segments (3->4)
	sta3temp_x = new vector<float>;
	sta3temp_y = new vector<float>;
	sta3temp_z = new vector<float>;
	sta4temp_x = new vector<float>;
	sta4temp_y = new vector<float>;
	sta4temp_z = new vector<float>;
	sta3dirtemp_x = new vector<float>;
	sta3dirtemp_y = new vector<float>;
	sta3dirtemp_z = new vector<float>;
	sta4dirtemp_x = new vector<float>;
	sta4dirtemp_y = new vector<float>;
	sta4dirtemp_z = new vector<float>;
	sizetemp = new vector<float>;
	int flag_34_1 = 0;
	for(int i4=sta4_x->size()-1 ; 0<=i4 ; i4--){ // Loop of all sta4 segments (2->4)
	  reVec34_x_34[i3][i4] = sta4_x->at(i4) - sta3_x->at(i3);
	  reVec34_y_34[i3][i4] = sta4_y->at(i4) - sta3_y->at(i3);
	  reVec34_z_34[i3][i4] = sta4_z->at(i4) - sta3_z->at(i3);
	  inner_product34_34[i3][i4] = (sta3dir_x->at(i3))*(reVec34_x_34[i3][i4])+(sta3dir_y->at(i3))*(reVec34_y_34[i3][i4])+(sta3dir_z->at(i3))*(reVec34_z_34[i3][i4]);
	  if( (inner_product34_34[i3][i4]>length(reVec34_x_34[i3][i4],reVec34_y_34[i3][i4],reVec34_z_34[i3][i4])*cos(60*3.14/180)) ){
	    cout<<"cone 3->4 (3->4)"<<endl;
	    flag_34_1 = 1;
	    sta3temp_x->push_back(sta3_x->at(i3));
	    sta3temp_y->push_back(sta3_y->at(i3));
	    sta3temp_z->push_back(sta3_z->at(i3));
	    sta4temp_x->push_back(sta4_x->at(i4));
	    sta4temp_y->push_back(sta4_y->at(i4));
	    sta4temp_z->push_back(sta4_z->at(i4));
	    sta3dirtemp_x->push_back(sta3dir_x->at(i3));
	    sta3dirtemp_y->push_back(sta3dir_y->at(i3));
	    sta3dirtemp_z->push_back(sta3dir_z->at(i3));
	    sta4dirtemp_x->push_back(sta4dir_x->at(i4));
	    sta4dirtemp_y->push_back(sta4dir_y->at(i4));
	    sta4dirtemp_z->push_back(sta4dir_z->at(i4));
	    sizetemp->push_back(1);
	    sizetemp->push_back(1);
	    sta4_x->erase(sta4_x->begin()+i4);
	    sta4_y->erase(sta4_y->begin()+i4);
	    sta4_z->erase(sta4_z->begin()+i4);
	    sta4dir_x->erase(sta4dir_x->begin()+i4);
	    sta4dir_y->erase(sta4dir_y->begin()+i4);
	    sta4dir_z->erase(sta4dir_z->begin()+i4);
	  }
	} // Loop of all sta4 segments (2->4)
	if(flag_34_1 == 1){
	  sta3_x->erase(sta3_x->begin()+i3);
	  sta3_y->erase(sta3_y->begin()+i3);
	  sta3_z->erase(sta3_z->begin()+i3);
	  sta3dir_x->erase(sta3dir_x->begin()+i3);
	  sta3dir_y->erase(sta3dir_y->begin()+i3);
	  sta3dir_z->erase(sta3dir_z->begin()+i3);
	}
	if(sizetemp->size()>1){
	  refFirstHit->push_back(muontemp_x->size());
	  if(sta3temp_x->size()>0){
	    for(unsigned int ii33=0 ; ii33<sta3temp_x->size() ; ii33++){
	      muontemp_x->push_back(sta3temp_x->at(ii33));
	      muontemp_y->push_back(sta3temp_y->at(ii33));
	      muontemp_z->push_back(sta3temp_z->at(ii33));
	      muondirtemp_x->push_back(sta3dirtemp_x->at(ii33));
	      muondirtemp_y->push_back(sta3dirtemp_y->at(ii33));
	      muondirtemp_z->push_back(sta3dirtemp_z->at(ii33));
	    }
	  }
	  if(sta4temp_x->size()>0){
	    for(unsigned int ii44=0 ; ii44<sta4temp_x->size() ; ii44++){
	      muontemp_x->push_back(sta4temp_x->at(ii44));
	      muontemp_y->push_back(sta4temp_y->at(ii44));
	      muontemp_z->push_back(sta4temp_z->at(ii44));
	      muondirtemp_x->push_back(sta4dirtemp_x->at(ii44));
	      muondirtemp_y->push_back(sta4dirtemp_y->at(ii44));
	      muondirtemp_z->push_back(sta4dirtemp_z->at(ii44));
	    }
	  }
	  nHits->push_back(sizetemp->size());
	}
	delete sta3temp_x;
	delete sta3temp_y;
	delete sta3temp_z;
	delete sta4temp_x;
	delete sta4temp_y;
	delete sta4temp_z;
	delete sta3dirtemp_x;
	delete sta3dirtemp_y;
	delete sta3dirtemp_z;
	delete sta4dirtemp_x;
	delete sta4dirtemp_y;
	delete sta4dirtemp_z;
	delete sizetemp;
      } // Loop of all sta3 segments (3->4)
    }
    ////////////////////// sta3->sta4 ////////////////////////////////////// 
    cout<<"---------------------------After 3->4-----------------------------------"<<endl;
    if(sta1_x->size()>0){
      cout<<"sta1_x = ( ";
      for(unsigned int sta1 = 0 ; sta1<sta1_x->size() ; sta1++){
	cout<<sta1_x->at(sta1)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta1_x->size() == 0){
      cout<<"sta1_x : No Segment"<<endl;
    }
    
    if(sta2_x->size()>0){
      cout<<"sta2_x = ( ";
      for(unsigned int sta2 = 0 ; sta2<sta2_x->size() ; sta2++){
	cout<<sta2_x->at(sta2)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta2_x->size() == 0){
      cout<<"sta2_x : No Segment"<<endl;
    }
    
    if(sta3_x->size()>0){
      cout<<"sta3_x = ( ";
      for(unsigned int sta3 = 0 ; sta3<sta3_x->size() ; sta3++){
	cout<<sta3_x->at(sta3)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta3_x->size() == 0){
      cout<<"sta3_x : No Segment"<<endl;
    }
    
    if(sta4_x->size()>0){
      cout<<"sta4_x = ( ";
      for(unsigned int sta4 = 0 ; sta4<sta4_x->size() ; sta4++){
	cout<<sta4_x->at(sta4)<<" ";
      }
      cout<<")"<<endl;
    } else if(sta4_x->size() == 0){
      cout<<"sta4_x : No Segment"<<endl;
    }
    cout<<"--------------------------------------------------------------"<<endl;
    
    eventsTree->Fill();
    
    delete sta1_x;
    delete sta1_y;
    delete sta1_z;
    delete sta2_x;
    delete sta2_y;
    delete sta2_z;
    delete sta3_x;
    delete sta3_y;
    delete sta3_z;
    delete sta4_x;
    delete sta4_y;
    delete sta4_z;
    delete sta1dir_x;
    delete sta1dir_y;
    delete sta1dir_z;
    delete sta2dir_x;
    delete sta2dir_y;
    delete sta2dir_z;
    delete sta3dir_x;
    delete sta3dir_y;
    delete sta3dir_z;
    delete sta4dir_x;
    delete sta4dir_y;
    delete sta4dir_z;
    delete muontemp_x;
    delete muontemp_y;
    delete muontemp_z;
    delete muondirtemp_x;
    delete muondirtemp_y;
    delete muondirtemp_z;
    delete refFirstHit;
    delete nHits;        
  } // Loop of all events
  myFile->cd();
  eventsTree->Write();
  myFile->Close();
}
