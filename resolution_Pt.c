void resolution_Pt()
{
	//TCanvas *c1 = new TCanvas();
	

	 gROOT->SetStyle("Plain");
     gStyle->SetLineWidth(3);
     //gStyle->SetOptStat("nem");
     gStyle->SetPadLeftMargin(0.14);

	TFile *input = new TFile("Q2=1_10x275_v5_new.root");	

	TTree *tree = (TTree*)input->Get("h7");

	Double_t pperp_gen,pperp_rec, eta_gen, eta_rec, tarberutyun;
	Int_t pid_gen, pid_rec;
//    double stdy[10]={0.0265,0.04448,0.05429,0.0525,0.04054,0.04249,0.05388,0.05502,0.04474,0.02673};
  //  double stdx[10]={0.1777,0.1772,0.1797,0.1786,0.1834,0.1836,0.1775,0.1796,0.177,0.1773};

	tree->SetBranchAddress("pperp_gen",&pperp_gen);
	tree->SetBranchAddress("pperp_rec",&pperp_rec);
    tree->SetBranchAddress("pid_gen",&pid_gen);
    tree->SetBranchAddress("eta_gen",&eta_gen);
    tree->SetBranchAddress("pid_rec",&pid_rec);
    tree->SetBranchAddress("eta_rec",&eta_rec);

    //TH1D *tarphi=new TH1D("tarphi", "", 500,-0.03,0.03 );
    //TH2D *tarphi=new TH2D("tarphi", "Resolution", 100,-3.14,3.14,100, -0.5,0.5 );
    //TGraph *gr= new TGraph(10,stdx,stdy);
    TH1D *tarphi[10];
    for(int i=0; i<10; i++){
        tarphi[i] = new TH1D(Form("tarphi_%d", i), Form("Resolution_%d", i), 100, -0.5, 0.5);
    }

	int entries = tree->GetEntries();
    int bin;


	for(int i=0;i<entries;i++){

		tree->GetEntry(i);
		if(pperp_rec<=4 && pperp_gen<=4){
            tarberutyun= pperp_rec - pperp_gen;
            //cout << tarberutyun << endl;
		  if(eta_rec > -3.5 && eta_rec < 3.5){
			 if(abs(pid_rec) == 211 /*|| abs(pid_rec)==321 */){
                    bin = int( (pperp_rec)/(0.4) );
                   // cout<<pperp_rec<<endl;
				    tarphi[bin]->Fill(tarberutyun);

		
			}
			
		  }
        }
		
	}

    double stdy[10], pperp_x[10], sxaly[10];
    double delta = 0.4;
    for(int i = 0; i < 10; i++){
        stdy[i] = tarphi[i]->GetRMS();
        sxaly[i] = tarphi[i]->GetRMSError();
        pperp_x[i] = (0 + delta/2 + delta*i);
        cout<<pperp_x[i]<<endl;

    }

    TGraphErrors *gr = new TGraphErrors(10, pperp_x, stdy, NULL, sxaly);


    TCanvas *c2 = new TCanvas("c2", "c2", 600, 450);
    c2->cd();
    gPad->SetBottomMargin(0.12);
    gr->SetMarkerSize(1.3);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(2);
    gr->GetYaxis()->SetTitle("Resolution #Delta P_{t}");
    gr->GetXaxis()->SetTitle("P_{t}");
    gr->GetYaxis()->SetTitleOffset(1);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetLimits(0,0.1);
    gr->GetYaxis()->SetRangeUser(0,0.1);
    //gr->GetYaxis()->SetTickLength(0);
    gr->GetYaxis()->SetNdivisions(80403);
    gr->SetTitle("");
    gr->Draw("AP");

    TCanvas *c1 = new TCanvas("c1","c1", 1800,900);
    c1->Divide(5,2);
    for(int i=0; i<10; i++){
        c1->cd(i+1);
        tarphi[i]->DrawCopy();
    }


    /*
    tarphi->SetMarkerSize(0.7);
    tarphi->SetMarkerStyle(20);
    tarphi->SetMarkerColor(4);
    tarphi->GetYaxis()->SetTitle("counts");
    tarphi->GetXaxis()->SetTitle("Phi_rec-Phi_gen");
    tarphi->GetYaxis()->SetTitleOffset(1.4);
    tarphi->Draw(); */
    /*gr->GetYaxis()->SetTitle("STD y");
    gr->GetXaxis()->SetTitle("STD x");
    gr->Draw("A*");*/


   /* c1->cd();
   pi0_q2_xb->SetLineColor(kBlue);
   pi0_q2_xb->SetLineWidth(2);
   pi0_q2_xb->GetXaxis()->SetTitle("x_{bj}");
   pi0_q2_xb->GetYaxis()->SetTitle("Q^{2}(Gev^{2})");
   pi0_q2_xb->GetYaxis()->SetTitleOffset(1.4);
   pi0_q2_xb->DrawClone("colz");
   c1->SaveAs("results/pi0_q2_xbj.pdf");

   c2->cd();
   pi0_pt_z->SetLineColor(kBlue);
   pi0_pt_z->SetLineWidth(2);
   pi0_pt_z->GetXaxis()->SetTitle("Z");
   pi0_pt_z->GetYaxis()->SetTitle("P_{T}");
   pi0_pt_z->GetYaxis()->SetTitleOffset(1.4);
   pi0_pt_z->DrawClone("colz");
   c2->SaveAs("results/pi0_pt_z.pdf");

   c3->cd();
   pi_q2_xb->SetLineColor(kBlue);
   pi_q2_xb->SetLineWidth(2);
   pi_q2_xb->GetXaxis()->SetTitle("x_{bj}");
   pi_q2_xb->GetYaxis()->SetTitle("Q^{2}(Gev^{2})");
   pi_q2_xb->GetYaxis()->SetTitleOffset(1.4);
   pi_q2_xb->DrawClone("colz");
   c3->SaveAs("results/pi+-_q2_xbj.pdf");

   c4->cd();
   pi_pt_z->SetLineColor(kBlue);
   pi_pt_z->SetLineWidth(2);
   pi_pt_z->GetXaxis()->SetTitle("Z");
   pi_pt_z->GetYaxis()->SetTitle("P_{T}");
   pi_pt_z->GetYaxis()->SetTitleOffset(1.4);
   pi_pt_z->DrawClone("colz");
   c4->SaveAs("results/pi+-_pt_z.pdf");

   c5->cd();
   jpsi_pt_z->SetLineColor(kBlue);
   jpsi_pt_z->SetLineWidth(2);
   jpsi_pt_z->GetXaxis()->SetTitle("Z");
   jpsi_pt_z->GetYaxis()->SetTitle("P_{T}");
   jpsi_pt_z->GetYaxis()->SetTitleOffset(1.4);
   jpsi_pt_z->DrawClone("colz");
   c5->SaveAs("results/jpsi_pt_z.pdf");

   c6->cd();
   jpsi_q2_xb->SetLineColor(kBlue);
   jpsi_q2_xb->SetLineWidth(2);
   jpsi_q2_xb->GetXaxis()->SetTitle("x_{bj}");
   jpsi_q2_xb->GetYaxis()->SetTitle("Q^{2}(Gev^{2})");
   jpsi_q2_xb->GetYaxis()->SetTitleOffset(1.4);
   jpsi_q2_xb->DrawClone("colz");
   c6->SaveAs("results/jpsi_q2_xb.pdf");
*/

}

