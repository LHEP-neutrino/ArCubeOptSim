
#include <iostream>
#include "TClassTable.h"

#include <vector>
#include <map>
#include <set>
#include <cmath>

#include <stdio>
//#include <stdint.h>
//#include <string>
//#include <stdlib>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"
#include "TROOT.h"


void pdesim_analysis() {
   //Get old file, old tree and set top branch address
   TFile *file = new TFile("./events.root");
   TTree *tr = (TTree*)file->Get("t1");
   Int_t neve = tr->GetEntries();
   
   Int_t eventid;
   Float_t prim_Xpos;
   Float_t prim_Ypos;
   Float_t prim_Zpos;
   Int_t totsteps;
   vector<Long64_t> trackid;
   vector<Long64_t> parentid;
   vector<Long64_t> firstparentid;
     
  tr->SetBranchAddress("eventid",&eventid); //ev number
  tr->SetBranchAddress("prim_Xpos",&prim_Xpos);
  tr->SetBranchAddress("prim_Ypos",&prim_Ypos);
  tr->SetBranchAddress("prim_Zpos",&prim_Zpos);
  tr->SetBranchAddress("totsteps",&totsteps);
  tr->SetBranchAddress("trackid",&trackid);
  tr->SetBranchAddress("parentid",&parentid);
  tr->SetBranchAddress("firstparentid",&firstparentid);
  
   for (Long64_t i=0;i<nentries; i++) {
      tr->GetEntry(i);
      if (run==runnr) newtree->Fill();
   }
   newtree->Print();
   newtree->AutoSave();
   delete oldfile;
   delete newfile;
}

