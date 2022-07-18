
void plot_data_220718() {

    TFile* ifile = new TFile("~/CERN/data/pulse_information-220-plus.root", "READ");
    // TTree* tree = (TTree*)ifile->Get("pulse_information");
    TH1D* hTOF_cut = (TH1D*)ifile->Get("TOF_CFD_veto");
    TH1D* hACT23_TOFc = (TH1D*)ifile->Get("hACT23_TOF");
    hTOF_cut->SetMarkerStyle(20);
    hTOF_cut->SetMarkerSize(1);
    hTOF_cut->SetMarkerColor(1);
    hTOF_cut->SetLineColor(1);
    hTOF_cut->SetTitle("+220 Mev, ACT01<50mV, ACT23<20mV, ACT45<14mV");
    // hTOF_cut->Draw("PE");

    // hACT23_TOFc->Draw("colz");

    // 220 MeV
    TF1* g1 = new TF1("1efit", "gaus", 38., 40.2);
    TF1* g2 = new TF1("2efit", "gaus", 40.2, 41.2);
    TF1* g3 = new TF1("3efit", "gaus", 41.2, 42);
    TF1* total = new TF1("mstotal","gaus(0)+gaus(3)+gaus(6)",38,42);

    double par[9];
    hTOF_cut->Fit(g1,"R 0");
    hTOF_cut->Fit(g2,"R+ 0");
    hTOF_cut->Fit(g3,"R+ 0");

    // Get the parameters from the fit
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);

    // Use the parameters on the sum
    total->SetParameters(par);
    hTOF_cut->Fit(total,"R+ 0");

    // Separate gaussian
    TF1* g1s = new TF1("1efits", "gaus", 38., 42);
    TF1* g2s = new TF1("2efits", "gaus", 38., 42);
    TF1* g3s = new TF1("3efits", "gaus", 38., 42);
    g1s->SetParameters(&par[0]);
    g2s->SetParameters(&par[3]);
    g3s->SetParameters(&par[6]);
    g1s->Draw("SAME");
    g2s->Draw("SAME");
    g3s->Draw("SAME");

    // TF1* total = new TF1("mstotal","gaus(0)+gaus(3)+gaus(6)",38,42);
    // hTOF_cut->Fit("efit");

}