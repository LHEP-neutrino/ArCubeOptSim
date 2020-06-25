
#include <iostream>
#include "TClassTable.h"

#include <vector>
#include <map>
#include <set>
#include <cmath>

//#include <stdio>
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


void pdesim_analysis(const char* filename,Double_t nrprim) {
   //Get old file, old tree and set top branch address
   TFile *file = new TFile(filename);
   TTree *tr = (TTree*)file->Get("t1");
   Int_t neve = tr->GetEntries();
   
   Int_t eventid;
   Double_t prim_Xpos;
   Double_t prim_Ypos;
   Double_t prim_Zpos;
   Int_t totsteps;
   vector<Long64_t> *trackid=0;
   vector<Long64_t> *firstparentid=0;
   vector<Int_t> evhits (neve,0);
     
  tr->SetBranchAddress("EvId",&eventid); //ev number
  tr->SetBranchAddress("prim_Xpos",&prim_Xpos);
  tr->SetBranchAddress("prim_Ypos",&prim_Ypos);
  tr->SetBranchAddress("prim_Zpos",&prim_Zpos);
  tr->SetBranchAddress("totsteps",&totsteps);
  tr->SetBranchAddress("trackid",&trackid);
  tr->SetBranchAddress("firstparentid",&firstparentid);

  Double_t pde_result = 0;
  
   for (Long64_t i=0;i<neve; i++) {
      tr->GetEntry(i);
      Long64_t evlen = firstparentid->size();
      Int_t temp_pid = 0;
      Int_t temp_pid_old = 0;
      Int_t nrhits = 0;
	 if(evlen!=0){
		 for (Long64_t j=0;j<evlen; j++) {
			 //std::cout << (*firstparentid)[j] << " ";
			 temp_pid = (*firstparentid)[j];
			 if((temp_pid!=temp_pid_old)){
				nrhits++;}// else std::cout <<temp_pid << "\n";
			 if(temp_pid==0) std::cout << "WARNING pid=0\n";
			 temp_pid_old = temp_pid;
	  	}
	 }
      //std::cout << "nrhits: " << nrhits << " totsteps: " << totsteps << std::endl;
      pde_result += nrhits/nrprim*0.38;
   }
   pde_result /= (Double_t)neve;
   pde_result *=100.0;
   std::cout << "Overall PDE: " << pde_result << "%" << std::endl;
   return;
}

