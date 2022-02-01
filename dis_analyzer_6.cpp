
#include <stdio.h>

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2D.h>
#include <TChain.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TStyle.h>
#include"stdio.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif


double vdot(double vv1[], double vv2[])
{
    double vdot;
    vdot=vv1[0]*vv2[0]+vv1[1]*vv2[1]+vv1[2]*vv2[2];
    return vdot;
}

double vmod(double vv[])
{
    double vmod;
    vmod=sqrt(vdot(vv,vv));
    return vmod;
}

void vprod(double vv1[],double vv2[],double vv3[])
{
    vv3[0]=vv1[1]*vv2[2]-vv1[2]*vv2[1];
    vv3[1]=vv1[2]*vv2[0]-vv1[0]*vv2[2];
    vv3[2]=vv1[0]*vv2[1]-vv1[1]*vv2[0];
    return;
}

int cos_check(double cosangle)
{
    int ok=1;
    if(fabs(cosangle)>1.0){
        if(fabs(cosangle)<1.00001){
            cosangle=cosangle/fabs(cosangle);
            ok=1;
        }
        else{
            ok=0;
        }
    }
    return ok;
}

double ang_betw_planes(double v1[], double w[], double v2[])
{
    double ang,cosang,wdotww;
    int ok1;
    double vnorm1[3],vnorm2[3],ww[3];
    ang=9999.0;
    vprod(v1,w,vnorm1);
    vprod(v2,w,vnorm2);
    cosang=(vdot(vnorm1,vnorm2))/(vmod(vnorm1)*vmod(vnorm2));
    ok1=cos_check(cosang);
    if(ok1==1){
        ang=acos(cosang);
    }
    vprod(vnorm1,vnorm2,ww);
    wdotww=vdot(w,ww);
    if(wdotww<0.0){
        ang=-ang;
    }
    return ang;
}

double ang_betw_planes_1(TVector3 v1, TVector3 w, TVector3 v2)
{
    double ang,cosang,wdotww;
    int ok1;
    TVector3 vnorm1;
    TVector3 vnorm2;
    TVector3 ww;
    ang=9999.0;
    vnorm1 = v1.Cross(w);
    vnorm2 = v2.Cross(w);
//    vprod(v1,w,vnorm1);
//    vprod(v2,w,vnorm2);
    cosang=(vnorm1.Dot(vnorm2)/(vnorm1.Mag() * vnorm2.Mag()) );
//    cosang=(vdot(vnorm1,vnorm2))/(vmod(vnorm1)*vmod(vnorm2));
    ok1=cos_check(cosang);
    if(ok1==1){
        ang=acos(cosang);
    }
    ww = vnorm1.Cross(vnorm2);
//    vprod(vnorm1,vnorm2,ww);
    wdotww=w.Dot(ww);
//    wdotww=vdot(w,ww);
    if(wdotww<0.0){
        ang=-ang;
    }
    return ang;
}


double get_Q2( TLorentzVector kb, TLorentzVector ks)
{
    return 2.*kb.Dot(ks);
}

double get_xbj( TLorentzVector kb, TLorentzVector ks, TLorentzVector pb)
{
    double q2 = get_Q2(kb, ks);
    return q2/(-2.0*(pb.Dot(ks)-pb.Dot(kb)));
}

double get_pt( TLorentzVector kb, TLorentzVector ks, TLorentzVector kh)
{
    TLorentzVector qprime = kb - ks;
    TVector3 qqprime = qprime.Vect();
    return kh.Perp(qqprime);
}

double get_z( TLorentzVector kb, TLorentzVector ks, TLorentzVector pb, TLorentzVector kh)
{
    TLorentzVector qprime = kb - ks;
    return ( pb.Dot(kh) )/( pb.Dot(qprime) );
}


double get_phi( TLorentzVector kb, TLorentzVector ks, TLorentzVector kh){
    TLorentzVector qprime = kb - ks;
    TVector3 qqprime = qprime.Vect();
    TVector3 vbeam = kb.Vect();
    TVector3 vhad = kh.Vect();
    return ang_betw_planes_1(vbeam,qqprime,vhad); 
}

double get_phis( TLorentzVector kb, TLorentzVector ks){
    TLorentzVector qprime = kb - ks;
    TVector3 qqprime = qprime.Vect();
    TVector3 vbeam = kb.Vect();
    TVector3 vhad(0,-1,0);
    return ang_betw_planes_1(vbeam,qqprime,vhad); 
}

double get_y( TLorentzVector kb, TLorentzVector ks, TLorentzVector pb)
{
    TLorentzVector qprime = kb - ks;
    return ( qprime.Dot(pb) )/( kb.Dot(pb) );
}

double get_s( TLorentzVector kb, TLorentzVector pb)
{
    TLorentzVector sum = kb + pb;
    return sum.Mag2();
}

double get_eps( TLorentzVector kb, TLorentzVector ks, TLorentzVector pb)
{
    double qq2 = get_Q2(kb, ks);
    double xbj = get_xbj(kb, ks, pb);
    double yy = get_y(kb, ks, pb);
    double mm = pb.Mag2();
    double gamma = 2.*mm*xbj/qq2;
    return (1.+yy - 0.25*gamma*gamma*yy*yy)/(1.-yy+0.25*yy*yy*(gamma*gamma+2.));
}

int get_bin_x( double xbj, double xmin, double k){
    int binx = 0; 
    double xtemp = xbj;
    while (xtemp/k > xmin){
        xtemp = xtemp/k;
        binx++ ; 
    }
    return binx;
}



int dis_analyzer_6(const char *fname="output18x275"){
    char inputfile[128];
    char outfile[128];
    sprintf(inputfile,"output/%s.root",fname);
    printf("%s\n",inputfile);
    gSystem->Load("libDelphes");
    TRandom3 ran3;

    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputfile);

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    
    // Get pointers to branches created by Delphes control card
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchEFlow = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    
    GenParticle * thisParticle; // Generated  particle branch
    Photon * thisPhoton;        //  Identified photon branch
    Electron * thisElectron;    //  Identified electron branch
    Jet * thisJet;
    // Track * thisEFlow; //?? Ask Miguel Arratia
    Track * thisTrack;
    GenParticle GParticle;



    // Initialize 4vectors and 3vectors of all particles
    TLorentzVector k4Beam, k4Scat, P4Beam, q4prime, p4Kaon, p4Pion, p4Prot, p4Had;
    TVector3 k3Beam, k3Scat, P3Beam, q3prime, p3Kaon, p3pion, p3Prot, p3Had;
    TLorentzVector k4BeamGen, k4ScatGen, P4BeamGen, q4primeGen, p4KaonGen, p4PionGen, p4ProtGen, p4HadGen;
    Double_t mele = 0.000511;
    Double_t mpion = 0.13957;
    Double_t mkaon = 0.49368;
    Double_t mprot = 0.93827;
    Double_t mhad = mpion;
    Double_t Ebeam, Eele, Epbeam, Ehad;

    Double_t Q2_gen, xbj_gen, pperp_gen, z_gen;
    Int_t ele_count;
    Double_t phi_gen, phis_gen;
    Int_t pid_gen;

    Double_t s_gen, y_gen, eps_gen;
 
    Double_t hmom_gen, eta_gen;
    Double_t wei = 0;
    Int_t pol = 0;

    Double_t E_gen;

    Double_t E_all,eta_all , theta_all;
/*
    sprintf(outfile,"%s_new10.root",fname);
    TFile *c = new TFile(outfile,"recreate");
    TTree *h7 = new TTree("h7","tt7");

    h7->Branch("Q2_gen",&Q2_gen,"Q2_gen/D"); 
    h7->Branch("xbj_gen",&xbj_gen,"xbj_gen/D"); 
    h7->Branch("pperp_gen",&pperp_gen,"pperp_gen/D"); 
    h7->Branch("z_gen",&z_gen,"z_gen/D"); 
    h7->Branch("phi_gen",&phi_gen,"phi_gen/D"); 
    h7->Branch("phis_gen",&phis_gen,"phis_gen/D"); 
    h7->Branch("pid_gen",&pid_gen,"pid_gen/I"); 
    h7->Branch("s_gen",&s_gen,"s_gen/D"); 
    h7->Branch("y_gen",&y_gen,"y_gen/D"); 
    h7->Branch("eps_gen",&eps_gen,"eps_gen/D"); 
    h7->Branch("hmom_gen",&hmom_gen,"hmom_gen/D"); 
    h7->Branch("eta_gen",&eta_gen,"eta_gen/D"); 
    h7->Branch("E_gen",&E_gen,"E_gen/D");
*/
    sprintf(outfile,"%s_new_all.root",fname);
    TFile *c2 = new TFile(outfile,"recreate");
    TTree *ha = new TTree("ha","tt7");

    ha->Branch("eta_all",&eta_all,"eta_all/D"); 
    ha->Branch("E_all",&E_all,"E_all/D");
    ha->Branch("theta_all",&theta_all,"theta_all/D");





//    h7->Branch("wei",&wei,"wei/D"); 
//    h7->Branch("pol",&pol,"pol/I"); 

//    k4Beam.SetPxPyPzE(0.,0.,18.,18);

    const int nbin_xA = 12;
    double xA_fact = pow(10.,-1./5.);
    double xA_BinLowEdge[nbin_xA+1];
    double xxA=0.1;
  //  printf("xBj low edge" );
    for (int ibin=0; ibin<=nbin_xA; ibin++) {
      //  printf("%8.3g, ", xxA );
        xA_BinLowEdge[nbin_xA-ibin] = xxA;
        xxA *= xA_fact;
    }
    //printf(" \n");
    const double Q2Min = 1.0, Q2Max=23.7137302;
    const int nbin_Q2 = 11;
    double Q2_BinLowEdge[nbin_Q2+1];
    double Q2_fact = pow(Q2Max/Q2Min,1.0/(double)nbin_Q2);
    double Q2val = Q2Min;
    for (int ibin=0; ibin<=nbin_Q2+1; ibin++) {
       // printf("%8.3g, ", Q2val );
        Q2_BinLowEdge[ibin] = Q2val;
        Q2val*=Q2_fact;
    }
    //printf("\n");
    
    int nbin_pt = 14;
    int nbin_z = 13;

    double pt_BinLowEdge[15]={0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 4.0};
    double z_BinLowEdge[14] = {0. ,0.05 ,0.1 ,0.15 ,0.2 ,0.25 ,0.3 ,0.4 ,0.5 ,0.6 ,0.7 ,0.8 ,0.9, 1.0};


/*
    for(int qbin = 1; qbin<nbin_Q2; ++qbin){
        for(int xbin = 1; xbin<nbin_xA; ++xbin){
            for(int zbin = 1; zbin<nbin_z; ++zbin){
                for(int ptbin = 1; ptbin<nbin_pt; ++ptbin){
                    printf("%9.5lf %9.5lf", Q2_BinLowEdge[qbin], Q2_BinLowEdge[qbin+1]);
                    printf("%9.5lf %9.5lf", xA_BinLowEdge[xbin], xA_BinLowEdge[xbin+1]);
                    printf("%6.2lf %6.2lf", z_BinLowEdge[zbin], z_BinLowEdge[zbin+1]);
                    printf("%6.2lf %6.2lf\n", pt_BinLowEdge[ptbin], pt_BinLowEdge[ptbin+1]);
                }
            }
        }
    }
*/
     

    TH2D *hh2[12][11];
    for(int ii = 0; ii<12 ; ++ii){
        for(int jj = 0; jj<11 ; ++jj){
            hh2[ii][jj] = new TH2D(Form("q2_%d_xbj_%d",jj,ii),Form("q2_%d_xbj_%d",jj,ii),nbin_z,z_BinLowEdge,nbin_pt,pt_BinLowEdge);
        }
    }


    TH2D *h2_q2_xb = new TH2D("q2_xb","q2_xb",nbin_xA,xA_BinLowEdge,nbin_Q2,Q2_BinLowEdge);

   // Loop over all events
    Int_t entry=0;
    Int_t nElectron=0;
    Int_t nJet=0;
    Int_t nPhoton=0;

    Int_t gen_pions, gen_kaons, gen_protons;

    Double_t rele_mom, gele_mom, ghad_mom;
    int bin, bin_h;
    TRef nn;
    TRef mm;

    int bin_xb, bin_q2;

    for(entry = 0; entry < numberOfEntries; ++entry)
//    for(entry = 0; entry < 100; ++entry)
    {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        gen_pions = 0;
        gen_kaons = 0;

        vector<TLorentzVector> gpion;

        vector<TLorentzVector> gkaon;

        vector<TLorentzVector> ghadron;
        vector<int> ghadron_pid;

        TLorentzVector all_p;
        vector<TLorentzVector> all_part;

        //  branchParticle holds the generated particles
        gele_mom = 0;
        for (int jpart=0; jpart<branchParticle->GetEntries(); jpart++){
            thisParticle = (GenParticle *) branchParticle->At(jpart);


            if (thisParticle->Status==4) {
                if (thisParticle->PID == 11) {
                    k4BeamGen = thisParticle->P4();
                }
                else if (thisParticle->PID == 2212) {
                    P4BeamGen = thisParticle->P4();
                }
            }
            if (thisParticle->Status==2) {
                if (abs(thisParticle->PID) == 111){
                    p4HadGen = thisParticle->P4();
                    ghadron.push_back(p4HadGen);
                    ghadron_pid.push_back(thisParticle->PID);
                }      
            }
            
            if(abs(thisParticle->PID) == 211){
                    p4HadGen = thisParticle->P4();
                    ghadron.push_back(p4HadGen);
                    ghadron_pid.push_back(thisParticle->PID);
                }

            if(thisParticle->PID == 11){
                //if(thisParticle->P > rele_momik){
                    all_p = thisParticle->P4();
                    all_part.push_back(all_p);
                //}
                        
            }
            

        }

        //  branchTrack holds the reconstructed particles
        rele_mom = 0;
        for (int iTrack = 0; iTrack < branchTrack->GetEntries(); ++iTrack ){
            thisTrack = (Track *) branchTrack->At(iTrack);
            if (thisTrack->PID == 11){
                if ( thisTrack->P > rele_mom){
                    k4Scat = thisTrack->P4();
                    rele_mom = thisTrack->P;
                    nn = thisTrack->Particle;
                    thisParticle = (GenParticle*) thisTrack->Particle.GetObject();
                    k4ScatGen = thisParticle->P4();
                    Q2_gen = get_Q2(k4BeamGen,k4ScatGen);
                    xbj_gen = get_xbj(k4BeamGen,k4ScatGen,P4BeamGen);
                }


            }
        }

       
        Q2_gen = get_Q2(k4BeamGen,k4ScatGen);
        xbj_gen = get_xbj(k4BeamGen,k4ScatGen,P4BeamGen);
        s_gen = get_s(k4BeamGen, P4BeamGen); 
        y_gen = get_y(k4BeamGen, k4ScatGen, P4BeamGen);
        eps_gen = get_eps(k4BeamGen, k4ScatGen, P4BeamGen);

        
        for(int kk = 0; kk < all_part.size();kk++){
            E_all = all_part[kk].E();
            eta_all = all_part[kk].Eta();
            theta_all = 2*TMath::ATan(TMath::Exp(-1.0*eta_all*TMath::RadToDeg()));
            if(eta_all > -4 && eta_all < 4){
                ha->Fill();
            }
        }
        
        /*
        if(xbj_gen > xA_BinLowEdge[0] && xbj_gen <= 0.1 && Q2_gen>1. && Q2_gen<Q2_BinLowEdge[11]){
            for( int vv = 0; vv < ghadron.size(); ++vv){
                
//                bin_xb = get_bin_x(xbj_rec, xA_BinLowEdge[0], 1./xA_fact);
//                bin_q2 = get_bin_x(Q2_rec, 1., Q2_fact);
                pperp_gen = get_pt(k4BeamGen, k4ScatGen, ghadron[vv]);
                z_gen = get_z(k4BeamGen, k4ScatGen, P4BeamGen,ghadron[vv]);
                phi_gen = get_phi(k4BeamGen, k4ScatGen, ghadron[vv]);
                phis_gen = get_phis(k4BeamGen, k4ScatGen);
                pid_gen = ghadron_pid[vv];
                hmom_gen = ghadron[vv].P();
                eta_gen = ghadron[vv].Eta();
                E_gen = ghadron[vv].E();

                h7->Fill();
                
//                hh2[bin_xb][bin_q2]->Fill(z_rec,pperp_rec);
            }


        }
        */


   }
/*
   h7->Write("c", TObject::kOverwrite);
   c->Close();
*/
   ha->Write("", TObject::kOverwrite);
   c2->Close();
/*
     h2_q2_xb->GetXaxis()->SetTitle("Q^{2}");
     h2_q2_xb->GetYaxis()->SetTitle("X{bj.}");
     h2_q2_xb->Draw("colz");
*/
    return entry;
}
