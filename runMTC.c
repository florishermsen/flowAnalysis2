/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//////////                                         //////////
//////////                runMTC.c                 //////////
//////////                                         //////////
//////////    Original Macro by F.A.W. Hermsen     //////////
//////////                                         //////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


// Define Variables

//Define Centrality class
Int_t CClass = 0;

//Define V1
Double_t V1 = 0.0109;

//Set Rapidity Limits and Bins
Double_t etaRange = 0.8;
Int_t etaBins = 6;

//Set sample size steps and multiplier
Int_t SSSteps = 10;
Int_t SSMultiplier = 200000;


#include <ctime>
#include <string>
#include "Riostream.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

#include "AliFlowTrackSimple.h"

using std::endl;
using std::cout;

Int_t cycleReport = 2000;

void runMTC() {
	cout<<endl;
	cout<<"Starting Toy Monte-Carlo for v1="<< V1 <<endl;
	cout<<endl;
	gSystem->Sleep(500);

	// Configure rapidity profile

	TF1 *etaDistribution = new TF1("fEtaDistribution","1+[0]*x^2",-etaRange,etaRange); // % dip around 0
	etaDistribution->SetNpx(50);

   etaDistribution->SetParameter(0,0.056); //10-30
   if(CClass==1) {
     	etaDistribution->SetParameter(0,0.061); //30-50
   } else if(CClass==2) {
     	etaDistribution->SetParameter(0,0.074); //60-80
   }

	// Configure fourrier transform of the azimuthal angle distribution
   TF1 *phiDistribution = new TF1("fPhiDistribution","1+2.*[0]*TMath::Cos(x)", 0, TMath::TwoPi());
   phiDistribution->SetNpx(50);
   phiDistribution->SetParName(0,"Directed Flow (v1)"); 
   phiDistribution->SetParameter(0,V1);

   // Configure fitting function for D+
   TF1 *v1FitPos = new TF1("v1FitPos","[0]*x", -etaRange, etaRange);
   v1FitPos->SetParameter(0,V1);

	// Configure fitting function for D-
   TF1 *v1FitNeg = new TF1("v1FitNeg","[0]*x", -etaRange, etaRange);
   v1FitNeg->SetParameter(0,-V1);

   // Create histogram for v1 vs N
   TGraphErrors *plotV1vsN = new TGraphErrors();
   
   // Create histogram for SoD vs N
   TGraph *plotSoDvsN = new TGraph();

	for(Int_t i=0; i<SSSteps; i++) {   

		//Sample Size
		Int_t sampleSize = (i+1)*SSMultiplier;

	   // Create histogram for D+
		TProfile *histFlowPos = new TProfile("histFlowPos","Directed Flow v_{1}(#eta) for D+",etaBins,-etaRange,etaRange);
	   histFlowPos->SetXTitle("#eta");
	   histFlowPos->SetYTitle("v_{1}");

		// Create histogram for D-
		TProfile *histFlowNeg = new TProfile("histFlowNeg","Directed Flow v_{1}(#eta) for D-",etaBins,-etaRange,etaRange);
	   histFlowNeg->SetXTitle("#eta");
	   histFlowNeg->SetYTitle("v_{1}");

		cout<<endl;
	   cout<<endl;
		cout<<"Starting with N="<< sampleSize <<endl;
		cout<<endl;
		gSystem->Sleep(500);

	   for(Int_t j=0; j<sampleSize; j++) {
	   	AliFlowTrackSimple *pTrack = new AliFlowTrackSimple();

			// Set Particle Charge
			pTrack->SetCharge((gRandom->Integer(2)>0.5 ? 1 : -1));

			// Set Particle PseudoRapidity
	   	pTrack->SetEta(etaDistribution->GetRandom());

	   	// Set v1 for this rapidity and get Azimutal Angle
	      phiDistribution->SetParameter(0,pTrack->Eta()*pTrack->Charge()*V1);
	      pTrack->SetPhi(phiDistribution->GetRandom());

	      // Fill the histrogram with v1 including RP Error
	      Double_t RPerror = 0.2;

			if(CClass==2) {
     			RPerror = 0.3; //60-80
   		}

	      if(pTrack->Charge() == 1){
	      	histFlowPos->Fill(pTrack->Eta(), TMath::Cos(pTrack->Phi()+gRandom->Gaus(0, RPerror*TMath::Pi())));
	      } else {
	      	histFlowNeg->Fill(pTrack->Eta(), TMath::Cos(pTrack->Phi()+gRandom->Gaus(0, RPerror*TMath::Pi())));
	      }

	      // Delete the track to prevent memory leak
	   	delete pTrack;

	   	// Output Cycle Report
		   if((j % cycleReport) == 0) {
		      cout <<"  .... "<< j << " particles processed ...."<<endl;
		   }
	   }

	   // retrieve a value for v1 from fit
	    // D+
	   cout<<endl;
	   histFlowPos->Fit("v1FitPos","Q");
	    // D-
	   cout<<endl;
	   histFlowNeg->Fit("v1FitNeg","Q");

	   // Get values for v1
	    // D+
	   Double_t v1Pos = v1FitPos->GetParameter(0);
	   Double_t v1PosError = v1FitPos->GetParError(0);
	    // D-
	   Double_t v1Neg = v1FitNeg->GetParameter(0);
	   Double_t v1NegError = v1FitNeg->GetParError(0);

	   // Calculate the mean v1 and its error
	   Double_t avgV1 = (v1Pos - v1Neg)/2;
	   Double_t avgV1Error = TMath::Sqrt(pow(v1PosError, 2)+pow(v1NegError, 2))/2;

	   // Calculate the SoD
	   Double_t SoD = avgV1/avgV1Error;

	   //Report Values
	   cout<<endl;
	   cout<<"Done!"<<endl;
	   cout<<endl;
	   cout<<"Number of D+: "<<histFlowPos->GetEntries()<<endl;
	   cout<<"v1 = "<< v1Pos <<" ± "<< v1PosError <<endl;
	   cout<<endl;
	   cout<<"Number of D-: "<<histFlowNeg->GetEntries()<<endl;
	   cout<<"v1 = "<< v1Neg <<" ± "<< v1NegError <<endl;
	   cout<<endl;
	   cout<<"Total v1 = "<< avgV1 <<" ± "<< avgV1Error <<endl;

	   plotV1vsN->SetPoint(i, sampleSize, avgV1);
	   plotV1vsN->SetPointError(i, 0, avgV1Error);

	   plotSoDvsN->SetPoint(i, sampleSize, SoD);

	   delete histFlowPos;
	   delete histFlowNeg;
	}

   TCanvas *c1 = new TCanvas("c1","v1 vs N",200,10,700,500);
   c1->SetGrid();
   plotV1vsN->Draw();

   cout<<endl;
   cout<<endl;
   cout<<endl;
	cout<<"v1 Evolution data:"<<endl;
	cout<<endl;
   plotV1vsN->Print("all");

   TCanvas *c2 = new TCanvas("c2","SoD vs N",200,10,700,500);
   c2->SetGrid();
   plotSoDvsN->Draw();


   cout<<endl;
   cout<<endl;
   cout<<endl;
	cout<<"SoD Evolution data:"<<endl;
	cout<<endl;
   plotSoDvsN->Print("all");
   cout<<endl;

   delete v1FitPos;
   delete v1FitNeg;

	delete phiDistribution;
   delete etaDistribution;
}