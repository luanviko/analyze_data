#include "stdlib.h"
#include <iostream>
#include "string.h"
#include <vector>
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"

// Function to find CFD time of pulse.
double CFD_timing(TH1D* waveform, double baseline, const int global_imin, const float startp, const float endp, const float percentage, double & CFD_amplitude) {
    double y_min = waveform->GetBinContent(global_imin)-baseline;
    double y_end = endp*y_min;
    double y_start = startp*y_min;
    double rise_amplitude = (endp-startp)*y_min;
    int j = global_imin;
    int j_end = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_end) && (j > 1)) {
        j--;
        j_end = j-1;
    }
    j = global_imin;
    int j_start = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_start) && (j > 1)) {
        j--;
        j_start = j+1;
    }
    double b = (y_end-y_start)/(double(j_end)-double(j_start));
    double a = y_start - j_start*b;
    double iCFD = (y_start+percentage*rise_amplitude-a)/b;
    CFD_amplitude = y_start+percentage*rise_amplitude;
    return iCFD;
}

// Function to find CFD time of pulse.
double CFD_timing_linearfit_wont_work(TH1D* waveform, double baseline, const int global_imin, const float startp, const float endp, const float percentage, double & CFD_amplitude) {
    double y_min = waveform->GetBinContent(global_imin)-baseline;
    double y_end = endp*y_min;
    double y_start = startp*y_min;
    double rise_amplitude = (endp-startp)*y_min;
    int j = global_imin;
    int j_end = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_end) && (j > 1)) {
        j--;
        j_end = j-1;
    }
    j = global_imin;
    int j_start = 0;
    while (((waveform->GetBinContent(j)-baseline) < y_start) && (j > 1)) {
        j--;
        j_start = j+1;
    }
    double b_int = (y_end-y_start)/(double(j_end)-double(j_start));
    double a_int = y_start - j_start*b_int;
    TF1* lfit = new TF1("lfit","pol1", j_start*0.625, (j_end+1)*0.625);
    lfit->SetParameters(-1000., 100);
    waveform->Fit("lfit");
    double a = lfit->GetParameter(0);
    double b = lfit->GetParameter(1);
    double iCFD_int = (y_start+percentage*rise_amplitude-a_int)/b_int;
    double iCFD = (y_start+percentage*rise_amplitude-a)/b;
    CFD_amplitude = y_start+percentage*rise_amplitude;
    std::cout << "a_int  a:"   << a_int << "  " << a  << std::endl;
    std::cout << "b_int  b:"   << b_int << "  " << b << std::endl;
    std::cout << "CFD_int CFD: " << iCFD_int << "  " << iCFD << std::endl;
    // std::cout << "j: " << j_end << "  " << j_start << std::endl;
    return iCFD;
}

// Fing time from global minimum
int global_timing(TH1D* waveform, double baseline, int s1, int s2) {
        double ymin = 99999999.;
        int imin = 0;
        int nsamples = waveform->GetNbinsX();
        for (int k = s1; k < s2; k++){
            if (waveform->GetBinContent(k)-baseline < ymin) {
                ymin = waveform->GetBinContent(k)-baseline;
                imin = k;
            }
        }
    return imin;
}

// Waveform baseline
double waveform_baseline(TH1D* waveform) {
    double sample_sum = 0.;
    int count = 0;
    for (int k=1; k <20; k++) {
        sample_sum += waveform->GetBinContent(k);
        count++;
    }
    double baseline = sample_sum/float(count);
    return baseline;
}

// Find charge of pulse in "width" before and after "center"
double pulse_charge(TH1D* waveform, double baseline, int center, int width) {
    double charge = 0.;
    for (int k = int(center-width); k <= int(center+width); k++)  charge += waveform->GetBinContent(k);
    return charge;
}


void save_all_waveforms2(TH1D* waveforms[], int global_timing[], double CFD_pos[], double CFD_amp[], int i, double ns_per_bin, double V_per_unit, double baselines[]) {
    
    // Convert units
    TCanvas* cc = new TCanvas(Form("cw%d",i),Form("Event %d",i),2000,1200);
    cc->Divide(4,4);
    for (int k=1; k<=16; k++) {
    for (int l=1; l<waveforms[k-1]->GetNbinsX()+1; l++) waveforms[k-1]->SetBinContent(l, (waveforms[k-1]->GetBinContent(l)-baselines[k-1])*V_per_unit );
    cc->SetTitle(Form("CH%d",k-1));
    cc->cd(k);

    // Plot global
    TGraph* graph1 = new TGraph();
    graph1->Set(1);
    graph1->SetMarkerSize(1);
    graph1->SetMarkerStyle(8);
    graph1->SetMarkerColor(2);
    graph1->SetPoint(0, global_timing[k-1]*ns_per_bin-0.5*ns_per_bin, waveforms[k-1]->GetBinContent(global_timing[k-1]));

    // Plot CFD
    TGraph* graph2 = new TGraph();
    graph2->Set(1);
    graph2->SetMarkerSize(1);
    graph2->SetMarkerStyle(8);
    graph2->SetMarkerColor(8);
    graph2->SetPoint(0, CFD_pos[k-1]*ns_per_bin, CFD_amp[k-1]/1000.);

    // Plot waveform
    waveforms[k-1]->SetTitle(Form("waveform%d; Time (ns); Amplitude (V)", k));
    waveforms[k-1]->Draw("HIST");

    // Plot minima
    graph1->Draw("SAME P");  
    graph2->Draw("SAME P");

    // std::cout << k << "  " << global_timing[k-1]*ns_per_bin << "  " << CFD_pos[k-1] << std::endl;
    // std::cout << k << "  " << waveforms[k-1]->GetBinContent(global_timing[k-1]) << "  " << CFD_amp[k-1]/1000. << std::endl;
    }
    cc->SaveAs(Form("./waveform_dump/CFD_all_waveforms_%d.png",i));
}

void analyze_data_220718() {

    // 240 MeV for runs 134, 135 and 137.
    std::vector<std::string> filenames;
    // int energy = 300;
    // std::string polarity= "plus";
    // filenames.push_back("~/CERN/data/root_run_000141_0000_clean.root");  
    // filenames.push_back("~/CERN/data/root_run_000165_0000_clean.root");  
    // filenames.push_back("~/CERN/data/root_run_000140_0000_clean.root");  

    // int energy = 240;
    // std::string polarity= "plus";
    // filenames.push_back("~/CERN/data/root_run_000134_0000_clean.root");  
    // filenames.push_back("~/CERN/data/root_run_000135_0000_clean.root");  
    // filenames.push_back("~/CERN/data/root_run_000137_0000_clean.root");  
    // filenames.push_back("~/CERN/data/root_run_000162_0000_clean.root");  

    // int energy = "220"
    int energy = 220;
    std::string polarity = "minus";
    filenames.push_back("~/CERN/data/root_run_000144_0000_clean.root");
    filenames.push_back("~/CERN/data/root_run_000152_0000_clean.root");


    // filenames.push_back("~/CERN/data/root_run_000135_0000_clean.root");
    // filenames.push_back("~/CERN/data/root_run_000137_0000_clean.root");

    // Open output file
    char oname[1024];
    sprintf(oname, "~/CERN/data/test_pulse_information-%d-%s.root",energy,polarity.c_str());
    TFile* ofile = new TFile(oname,"RECREATE");
    TTree* otree = new TTree("pulse_information", "Pulse information");

    // Create output branches
    double oglobal_times[16];
    double oCFD_times[16];
    double obaselines[16];
    double oamplitudes[16];
    double ocharges[16];
    double oCFD_amplitudes[16];

    // Set addresses to new branch
    otree->Branch("global_times", &oglobal_times, "global_times/D");
    otree->Branch("CFD_times", &oCFD_times, "CFD_times/D");
    otree->Branch("amplitudes", &oamplitudes, "baseline/D");
    otree->Branch("charges", &ocharges, "charges/D");
    otree->Branch("CFD_amplitudes", &oCFD_amplitudes, "baseline/D");

    // Horizontal scale conversion factor
    double ns_per_bin = 0.625;
    double V_per_unit = 1./4096.*2.5*1.e3;

    int TOF_bins = 42;
    // int TOF_bins = 32;
    // int TOF_bins = 28;
    double TOF_min  = 38.;
    double TOF_max  = 42.;

    // Histograms to be filled in
    TH1D* hTOF_CFD    = new TH1D("TOF_CFD", "TOF; Time difference (ns); No. of events", TOF_bins, TOF_min, TOF_max);
    TH1D* hTOF_global = new TH1D("TOF_global", "TOF; Time difference (ns); No. of events", TOF_bins, TOF_min, TOF_max);
    TH1D* hTOF_cut    = new TH1D("TOF_CFD_veto", "TOF/ ACT01<50mV/ ACT23<20mV/ ACT45<14mV; Time difference (ns); No. of events", TOF_bins, TOF_min, TOF_max);
    TH2D* hACT01_LGC  = new TH2D("hACT01_LGC", "ACT01 by LGC; LGC amplitudes (mV); ACT01 amplitudes (mV)", 500, 0., 200., 500, 0., 3800.);
    TH2D* hACT01_TOF  = new TH2D("hACT01_TOF", "ACT01 by TOF; TOF (ns); ACT01 amplitudes (mV); TOF (ns)", TOF_bins, TOF_min, TOF_max, 500, 0., 3800.);
    TH2D* hACT23_TOF  = new TH2D("hACT23_TOF", "ACT23 by TOF; TOF (ns); ACT23 amplitudes (mV)", TOF_bins, TOF_min, TOF_max, 25, TOF_min, TOF_max);
    TH2D* hACT45_TOF  = new TH2D("hACT45_TOF", "ACT45 by TOF; TOF (ns); ACT45 amplitudes (mV)", TOF_bins, TOF_min, TOF_max, 500, 0., 1200.);
    TH2D* hACT01_TOFc = new TH2D("hACT01_TOFc", "ACT01 by TOF/ ACT01 veto; TOF (ns); ACT01 amplitudes (mV)", TOF_bins, TOF_min, TOF_max, 500, 0., 200.);
    TH2D* hACT23_TOFc = new TH2D("hACT23_TOFc", "ACT23 by TOF/ ACT01 veto; TOF (ns); ACT23 amplitudes (mV)", TOF_bins, TOF_min, TOF_max, 500, 0., 600.);
    TH2D* hACT45_TOFc = new TH2D("hACT45_TOFc", "ACT45 by TOF/ ACT01 veto; TOF (ns); ACT23 amplitudes (mV)", TOF_bins, TOF_min, TOF_max, 500, 0., 1200.);

    // To hold information from root file
    TH1D *waveforms[16];
    int freqsetting;

    //              0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15
    int s1[16] = {180, 160, 150, 150, 190, 170, 190, 150, 140, 140, 140, 140, 180, 180, 180, 180}; // ns
    int s2[16] = {280, 280, 250, 250, 290, 290, 280, 250, 220, 220, 220, 220, 250, 250, 250, 250}; // ns
    for (int i=0; i < 16; i++) s1[i]=s1[i]/ns_per_bin; // ns -> sample number. implicit int.
    for (int i=0; i < 16; i++) s2[i]=s2[i]/ns_per_bin; // ns -> sample number. implicit int.

    // Open trees
    TChain *digi1 = new TChain("midas_data1");
    TChain *digi2 = new TChain("midas_data2");  
    for(int i=0; i<(int)filenames.size(); i++){
        digi1->AddFile(filenames[i].c_str());
        digi2->AddFile(filenames[i].c_str());
    }  

    // Set branch addresses
    digi1->SetBranchAddress("freqsetting",&freqsetting);
    for(int i=0; i<8; i++){
        waveforms[i] = NULL;
        digi1->SetBranchAddress(Form("Channel%d",i),&(waveforms[i]));
    }  
    for(int i=0; i<8; i++){
        waveforms[i+8] = NULL;
        digi2->SetBranchAddress(Form("Channel%d",i),&(waveforms[i+8]));
    }  

    // Number of events
    Long64_t N;
    Long64_t N1 = digi1->GetEntries();
    Long64_t N2 = digi2->GetEntries();
    if (N1>=N2) N=N2;
    if (N1<N2) N=N1;

    // Do analysis
    std::cout << "Events to analyze: " << N << std::endl;
    for(Long64_t i=0; i<N; i++){

        // Show progress on terminal
        printf("Progress: %4.2f%%.\r", (double(i)+1)/double(N)*100.);
        
        digi1->GetEntry(i);
        digi2->GetEntry(i);

        int global_times[16];
        double CFD_times[16];
        double baselines[16];
        double amplitudes[16];
        double charges[16];
        double CFD_amplitudes[16];
        for(int j=0; j<16; j++){
            // Find values as arbitrary units
            baselines[j]    = waveform_baseline(waveforms[j]);
            global_times[j] = global_timing(waveforms[j], baselines[j], s1[j], s2[j]);
            CFD_times[j]    = CFD_timing(waveforms[j], baselines[j], global_times[j], 0.1, 0.9, 0.2, CFD_amplitudes[j]);
            charges[j]      = -1.*pulse_charge(waveforms[j], baselines[j], global_times[j], 20);
            amplitudes[j]   = -1.*(waveforms[j]->GetBinContent(global_times[j])-baselines[j]);
            // std::cout << amplitudes[j] << std::endl;
        }

        // // Uncomment if you want to save waveforms
        // // Waveforms already converted to ns and V.
        // save_all_waveforms2(waveforms, global_times, CFD_times, CFD_amplitudes,  i, ns_per_bin, V_per_unit, baselines);
    
        for (int j=0; j<16; j++) {
            // Convert units 
            baselines[j]      = baselines[j]*V_per_unit;
            global_times[j]   = global_times[j]*ns_per_bin;
            CFD_times[j]      = CFD_times[j]*ns_per_bin;
            charges[j]        = charges[j]*ns_per_bin*V_per_unit;
            amplitudes[j]     = amplitudes[j]*V_per_unit;
            CFD_amplitudes[j] = CFD_amplitudes[j]/1000.;

            // Fill branches
            obaselines[j]    = baselines[j];
            oglobal_times[j] = global_times[j]; 
            oCFD_times[j]    = CFD_times[j];
            ocharges[j]      = charges[j]; 
            oamplitudes[j]   = amplitudes[j];
            oCFD_amplitudes[j] = CFD_amplitudes[j];

            
        }

        otree->Fill();
    
        // Organizing data into the variables I want
        double TOF_CFD = 0.25*(-CFD_times[8]-CFD_times[9]-CFD_times[10]-CFD_times[11]+CFD_times[12]+CFD_times[13]+CFD_times[14]+CFD_times[15]);
        double TOF_global = 0.25*(-global_times[8]-global_times[9]-global_times[10]-global_times[11]+global_times[12]+global_times[13]+global_times[14]+global_times[15]); 
        double ACT01_charge = charges[0]+charges[1];
        double ACT23_charge = charges[2]+charges[3];
        double ACT45_charge = charges[4]+charges[5]; 
        double LGC_charge    = charges[6];
        double ACT01_amp = amplitudes[0]+amplitudes[1];
        double ACT23_amp = amplitudes[2]+amplitudes[3];
        double ACT45_amp = amplitudes[4]+amplitudes[5]; 
        double LGC_amp    = amplitudes[6];

        // Apply cut of <140 mV for ACT01.
        // Apply cut of <150 mV for ACT23.
        if (ACT01_amp < 50. && ACT23_amp < 20. && ACT45_amp < 14.) hTOF_cut->Fill(TOF_CFD);
        if (ACT01_amp < 50. && ACT23_amp < 20.) hACT01_TOFc->Fill(TOF_CFD, ACT01_amp);
        if (ACT01_amp < 50. && ACT23_amp < 20.) hACT23_TOFc->Fill(TOF_CFD, ACT23_amp);
        if (ACT01_amp < 50. && ACT23_amp < 20.) hACT45_TOFc->Fill(TOF_CFD, ACT45_amp);

        // Filling histograms with no cuts
        hTOF_CFD->Fill(TOF_CFD);
        hTOF_global->Fill(TOF_global);
        hACT01_LGC->Fill(LGC_amp,ACT01_amp);
        hACT01_TOF->Fill(TOF_CFD, ACT01_amp);
        hACT23_TOF->Fill(TOF_CFD, ACT23_amp);
        hACT45_TOF->Fill(TOF_CFD, ACT45_amp);

        // std:cout << ACT_01_amp << "  " << LGC_amp << std::endl;

        // std::cout << "TOF: " << TOF_CFD << std::endl;
        // std::cout << "Global: " << TOF_global << std::endl;
        // std::cout << "ACT_01: " << ACT_01_amp << std::endl;
        // std::cout << "ACT_23: " << ACT_01_amp << std::endl;
        // std::cout << "ACT_45: " << ACT_01_amp << std::endl;
    }

    ofile->Write();
    printf("Progress: 100.00%%\n");

    char ographs[2048];
    TCanvas* c1 = new TCanvas("c1", "TOF", 600, 500);
    c1->cd(0);
    hTOF_CFD->SetLineColor(2);
    hTOF_CFD->Draw("HIST");
    // hTOF_global->Draw("SAME");
    sprintf(ographs, "./images/%dMev-%s-hTOF_CFD.pdf", energy, polarity.c_str());
    c1->SaveAs(ographs);

    TCanvas*c2 = new TCanvas("c2", "TOF with cut", 700, 500);
    c2->cd(0);
    hTOF_cut->SetMarkerStyle(20);
    hTOF_cut->SetMarkerSize(1);
    hTOF_cut->SetMarkerColor(1);
    hTOF_cut->SetLineColor(1);
    hTOF_cut->Draw("PE0");
    
    // hTOF_CFD->Draw("SAME");
    sprintf(ographs, "./images/%dMev-%s-hTOF_cut.pdf", energy, polarity.c_str());
    c2->SaveAs(ographs);

    TCanvas*c3 = new TCanvas("c3", "", 600, 500);
    c3->cd(0);
    hACT01_TOF->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT01_TOF.pdf", energy, polarity.c_str());
    c3->SaveAs(ographs);

    TCanvas*c4 = new TCanvas("c4", "", 600, 500);
    c4->cd(0);
    hACT01_LGC->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT01_LGC.pdf", energy, polarity.c_str());
    c4->SaveAs(ographs);

    TCanvas*c5 = new TCanvas("c5", "", 600, 500);
    c5->cd(0);
    hACT23_TOF->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT23_TOF.pdf", energy, polarity.c_str());
    c5->SaveAs(ographs);

    TCanvas*c6 = new TCanvas("c6", "", 600, 500);
    c6->cd(0);
    hACT45_TOF->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT45_TOF.pdf", energy, polarity.c_str());
    c6->SaveAs(ographs);

    TCanvas*c7 = new TCanvas("c7", "", 600, 500);
    c7->cd(0);
    hACT01_TOFc->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT01_TOFc.pdf", energy, polarity.c_str());
    c7->SaveAs(ographs);

    TCanvas*c8 = new TCanvas("c8", "", 600, 500);
    c8->cd(0);
    hACT23_TOFc->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT23_TOFc.pdf", energy, polarity.c_str());
    c8->SaveAs(ographs);

    TCanvas*c9 = new TCanvas("c9", "", 600, 500);
    c9->cd(0);
    hACT45_TOFc->Draw("colz");
    sprintf(ographs, "./images/%dMev-%s-hACT45_TOFc.pdf", energy, polarity.c_str());
    c9->SaveAs(ographs);

    ofile->Close();
}