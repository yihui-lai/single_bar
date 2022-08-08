#include "TStyle.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include <TMath.h>
//#include <String.h>
#include <iostream>
#include <utility>

// UV and RGB SiPM QE
double UV_sipm_QE_x[23] = {361.161, 364.766, 379.794, 387.614,
                           396.624, 406.226, 411.617, 426.594, 436.769, 455.931, 477.492,
                           496.061, 517.627, 547.583, 573.349, 598.521, 615.299, 649.46,
                           671.039, 705.202, 755.548, 773.531, 798.108};
double UV_sipm_QE_y[23] = {0.770120854, 0.787348933, 0.879304547,
                           0.942520324, 0.982752141, 1, 0.982752141, 0.942520324, 0.890796527,
                           0.816088771, 0.741381015, 0.683901339, 0.620685563, 0.545977807,
                           0.488498131, 0.448266313, 0.413790375, 0.356330478, 0.32759064,
                           0.275866843, 0.201139308, 0.178155349, 0.149415511};
double RGB_sipm_QE_x[29] = {305.28, 318.47, 334.67, 352.06,
                            370.06, 396.44, 416.23, 443.81, 466.6, 477.39, 491.78, 515.17,
                            529.56, 556.53, 582.91, 610.49, 636.26, 663.24, 684.22, 712.39,
                            738.76, 755.55, 774.73, 795.11, 825.68, 850.26, 874.23, 894.61, 900.61};
double RGB_sipm_QE_y[29] = {0.034678173, 0.144499016, 0.271678829,
                            0.427750492, 0.525998688, 0.635839415, 0.705195761, 0.786124754,
                            0.87860651, 0.907518244, 0.936410093, 0.994213676, 1, 0.97687459,
                            0.942196417, 0.90173192, 0.849714661, 0.78033843, 0.734107494, 0.664731264,
                            0.583802271, 0.520232248, 0.485554075, 0.427750492, 0.364160585, 0.289017916,
                            0.225428009, 0.167624426, 0.144499016};
// make RGB, UV QE a function
double RGB_QE(double *x, double *par)
{
    int ntry = 28;
    if (x[0] >= RGB_sipm_QE_x[0] && x[0] < RGB_sipm_QE_x[ntry])
    {
        for (int i = 0; i < ntry; i++)
        {
            if (x[0] >= RGB_sipm_QE_x[i] && x[0] < RGB_sipm_QE_x[i + 1])
                return 0.325 * ((RGB_sipm_QE_y[i + 1] - RGB_sipm_QE_y[i]) * (x[0] - RGB_sipm_QE_x[i]) / (RGB_sipm_QE_x[i + 1] - RGB_sipm_QE_x[i]) + RGB_sipm_QE_y[i]);
        }
        return 0;
    }
    else
    {
        return 0;
    }
}
double UV_QE(double *x, double *par)
{
    int ntry = 22;
    if (x[0] >= UV_sipm_QE_x[0] && x[0] < UV_sipm_QE_x[ntry])
    {
        for (int i = 0; i < ntry; i++)
        {
            if (x[0] >= UV_sipm_QE_x[i] && x[0] < UV_sipm_QE_x[i + 1])
                return 0.43 * ((UV_sipm_QE_y[i + 1] - UV_sipm_QE_y[i]) * (x[0] - UV_sipm_QE_x[i]) / (UV_sipm_QE_x[i + 1] - UV_sipm_QE_x[i]) + UV_sipm_QE_y[i]);
        }
        return 0;
    }
    else
    {
        return 0;
    }
}


//BGO absorption, cm-1
double ABS_BGO_x[53] = {313.576, 313.931, 313.943,
                        313.958, 313.971, 314.294, 314.311, 314.326, 315.06, 315.075, 315.088, 315.45,
                        316.222, 316.605, 316.646, 317.797, 318.584, 320.176, 320.558, 322.147, 323.343,
                        324.943, 327.356, 329.765, 331.777, 335.401, 338.223, 341.046, 343.063, 345.079,
                        349.116, 352.75, 356.384, 360.424, 368.102, 377.801, 390.333, 403.674, 417.824,
                        433.997, 442.891, 453.402, 467.552, 483.724, 500.705, 515.26, 524.559, 540.327,
                        549.626, 564.181, 573.883, 585.206, 594.505};
double ABS_BGO_y[53] = {0.498203, 0.448958, 0.460516,
                        0.476094, 0.488657, 0.408256, 0.425341, 0.439411, 0.366047, 0.381123, 0.393685, 0.351476, 0.3158,
                        0.294193, 0.334897, 0.273594, 0.251989, 0.227371, 0.205262, 0.177629, 0.16105, 0.14397, 0.13041,
                        0.114338, 0.104797, 0.089733, 0.082205, 0.074677, 0.070161, 0.06464, 0.058624, 0.054113, 0.0491,
                        0.046601, 0.042104, 0.038117, 0.036149, 0.035188, 0.034231, 0.034285, 0.033812, 0.032842, 0.031382,
                        0.031436, 0.030488, 0.030537, 0.030568, 0.030621, 0.03015, 0.030198, 0.028723, 0.030269, 0.0303};

//PWO absorption, cm-1
double ABS_PWO_y[35] = {0.01, 0.02, 2.01,
                        4, 4.023302971, 4.255916788, 4.523924776, 4.815007415, 5.160917405, 5.560374547, 6.651855202,
                        7.330625889, 8.593280055, 9.708078092, 11.15523627, 13.0145633, 15.48275221, 18.70837387, 21.38077014, 23.02555837, 23.63172323, 23.9469336,
                        23.9469336, 24.27007742, 24.94449849, 24.94449849, 25.29596276, 25.29596276, 25.65681445, 25.65681445, 25.65681445, 25.65681445, 25.65681445, 26.41170567, 26.41170567};
double ABS_PWO_x[35] = {3.5, 3.44401, 3.41084, 3.37832, 3.376604418,
                        3.370894244, 3.365833707, 3.358266778, 3.353244055, 3.345733565, 3.330813064, 3.320945706, 3.303804576, 3.2868482, 3.270064988,
                        3.246373208, 3.220719845, 3.186369517, 3.13949987, 3.089725278, 3.037376425, 2.986779095, 2.937832906, 2.892332871, 2.846415974,
                        2.801927859, 2.758815138, 2.717003092, 2.678040871, 2.638629504, 2.600355854, 2.563181955, 2.52705078, 2.491924071, 1};

// make BGO, PWO absorption a function
double ABS_BGO(double *x, double *par)
{
    int ntry = 52;
    if (x[0] >= ABS_BGO_x[0] && x[0] < ABS_BGO_x[ntry])
    {
        for (int i = 0; i < ntry; i++)
        {
            if (x[0] >= ABS_BGO_x[i] && x[0] < ABS_BGO_x[i + 1])
                return 5 * ((ABS_BGO_y[i + 1] - ABS_BGO_y[i]) * (x[0] - ABS_BGO_x[i]) / (ABS_BGO_x[i + 1] - ABS_BGO_x[i]) + ABS_BGO_y[i]);
        }
        return 0;
    }
    else
    {
        return 0;
    }
}
double ABS_PWO(double *x, double *par)
{
    int ntry = 34;
    if (x[0] >= 1239.84187 / ABS_PWO_x[0] && x[0] < 1239.84187 / ABS_PWO_x[ntry])
    {
        for (int i = 0; i < ntry; i++)
        {
            if (x[0] >= 1239.84187 / ABS_PWO_x[i] && x[0] < 1239.84187 / ABS_PWO_x[i + 1])
                return 5 * ((1 / ABS_PWO_y[i + 1] - 1 / ABS_PWO_y[i]) * (x[0] - 1239.84187 / ABS_PWO_x[i]) / (1239.84187 / ABS_PWO_x[i + 1] - 1239.84187 / ABS_PWO_x[i]) + 1 / ABS_PWO_y[i]);
        }
        return 0;
    }
    else
    {
        return 0;
    }
}


double Get_Eff(TH1F *h_PhoSpectrum, bool use_UV, double wvmin, double wvmax)
{
    //get the detection effeciency in window (wvmin, wvmax)
    assert(wvmin <= wvmax);
    double wc_low = wvmin;
    double wc_high = wvmax;
    TF1 *rat_RGB_sipm_QE = new TF1("rat_RGB_sipm_QE", RGB_QE, wc_low, wc_high, 0);
    TF1 *rat_UV_sipm_QE = new TF1("rat_UV_sipm_QE", UV_QE, wc_low, wc_high, 0);
    double rat = 0;
    for (int ix = wc_low; ix < wc_high; ix++)
    {
        if (use_UV)
        {
            rat += rat_UV_sipm_QE->Eval(h_PhoSpectrum->GetBinCenter(ix)) * h_PhoSpectrum->GetBinContent(ix);
        }
        else
        {
            rat += rat_RGB_sipm_QE->Eval(h_PhoSpectrum->GetBinCenter(ix)) * h_PhoSpectrum->GetBinContent(ix);
        }
    }
    return rat / h_PhoSpectrum->Integral(300, 1000);

}

void Get_Nphoton(char *filename, double *output_, bool is_bgo, bool makeplots)
{
    //Get the number of photons

    TCanvas *c1 = new TCanvas();

    TFile *f = new TFile(filename);
    double Nphoton_generate[2]; // [S, C]
    TH1F *h_S_generate = new TH1F("h_S_generate", "", 10000, 0, 1000000);
    ((TTree *)f->Get("tree"))->Draw("ECAL_f_total_S>>h_S_generate");
    TH1F *h_C_generate = new TH1F("h_C_generate", "", 10000, 0, 200000);
    ((TTree *)f->Get("tree"))->Draw("ECAL_f_total_C>>h_C_generate");
    Nphoton_generate[0] = h_S_generate->GetMean();
    Nphoton_generate[1] = h_C_generate->GetMean();
    cout << "On average " << Nphoton_generate[0]<<" S, "<<Nphoton_generate[1]<<" C" << endl;


    double Nphoton_detect[4]; 
    TH1F *h_S_DetectFront = new TH1F("h_S_DetectFront", "", 5000, 0, Nphoton_generate[0] / 2); // detected photon will be much less than generated photons
    ((TTree *)f->Get("tree"))->Draw("SDdetected_ff_S>>h_S_DetectFront");
    cout << "S hit front: " << h_S_DetectFront->GetMean() << endl;
    TH1F *h_S_DetectRear = new TH1F("h_S_DetectRear", "", 5000, 0, Nphoton_generate[0] / 2);
    ((TTree *)f->Get("tree"))->Draw("SDdetected_rr_S>>h_S_DetectRear");
    cout << "S hit back: " << h_S_DetectRear->GetMean() << endl;
    TH1F *h_C_DetectFront = new TH1F("h_C_DetectFront", "", 5000, 0, Nphoton_generate[1] / 2);
    ((TTree *)f->Get("tree"))->Draw("SDdetected_ff_C>>h_C_DetectFront");
    cout << "C hit front: " << h_C_DetectFront->GetMean() << endl;
    TH1F *h_C_DetectRear = new TH1F("h_C_DetectRear", "", 5000, 0, Nphoton_generate[1] / 2);
    ((TTree *)f->Get("tree"))->Draw("SDdetected_rr_C>>h_C_DetectRear");
    cout << "C hit back: " << h_C_DetectRear->GetMean() << endl;
    Nphoton_detect[0] = h_S_DetectFront->GetMean();
    Nphoton_detect[1] = h_S_DetectRear->GetMean();
    Nphoton_detect[2] = h_C_DetectFront->GetMean();
    Nphoton_detect[3] = h_C_DetectRear->GetMean();

    // Use the spectrum (front and back should be the same) of collected C photons to calculate the SiPM eff.
    // Assume filter window is 550 nm to 1000nm
    double Efficiency =  Get_Eff((TH1F *)f->Get("h_phot_lambda_ECAL_f_collect_Ceren"), false, 550, 1000); 

    if(makeplots){
        
        TCanvas *c = new TCanvas();
        //c->SetLogy();
        gStyle->SetOptStat(0000);
        gStyle->SetOptTitle(0);
        TList *l = c->GetListOfPrimitives();
        TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
        p1->SetGrid();
        p1->SetLeftMargin(0.14);
        p1->SetRightMargin(0.14);
        p1->Draw();
        p1->cd();
    
    
        int num = ((TTree *)f->Get("tree"))->GetEntries();
        // Get the spectrum of generated S and C, normalize
        TH1F *hCeren = (TH1F *)f->Get("h_phot_lambda_ECAL_f_produce_Ceren");
        hCeren->Scale(1.0 / num);
        TH1F *hScinti = (TH1F *)f->Get("h_phot_lambda_ECAL_f_produce_Scin");
        hScinti->Scale(1.0 / num);
        double NS_max = hScinti->GetMaximum();
        double NC_max = hCeren->GetMaximum();
        hCeren->SetLineColor(kRed);
        hScinti->SetLineColor(kBlue);
        hCeren->SetLineWidth(3);
        hScinti->SetLineWidth(3);
        hScinti->GetXaxis()->SetTitle("Wavelength (nm)");
        hScinti->GetYaxis()->SetTitle("Number of photons (normalized)");
        hScinti->Draw("HIST");
        hCeren->Draw("same HIST");
    
        // Get the spectrum of S and C tha hits the sipm, normalize
        TH1F *hCeren_d = (TH1F *)f->Get("h_phot_lambda_ECAL_f_collect_Ceren");
        hCeren_d->Scale(1.0 / num);
        TH1F *hScinti_d = (TH1F *)f->Get("h_phot_lambda_ECAL_f_collect_Scin");
        hScinti_d->Scale(1.0 / num);
        double NS_det_max = hScinti_d->GetMaximum();
        double NC_det_max = hCeren_d->GetMaximum();
    
        hCeren_d->SetLineColor(kBlack);
        hScinti_d->SetLineColor(kOrange);
        hCeren_d->SetLineWidth(3);
        hScinti_d->SetLineWidth(3);
        hScinti_d->Draw("same HIST");
        hCeren_d->Draw("same HIST");
    
        TString text;
        text = Form("%5.1g", Nphoton_generate[0]);
        text.ReplaceAll("e+0","x10^{");
        text.Append("}");
        hScinti->SetTitle("Scintillation photons ( N_{Generated}= "+ text+")");
        text = Form("%5.1g", Nphoton_generate[1]);
        text.ReplaceAll("e+0","x10^{");
        text.Append("}");
        hCeren->SetTitle("Cerenkov photons ( N_{Generated}= "+ text+")");
        text = Form("%5.1g", Nphoton_detect[0]);
        text.ReplaceAll("e+0","x10^{");
        text.Append("}");
        hScinti_d->SetTitle("Scintillation photons ( N_{Reach_end}= "+ text+")");
        text = Form("%5.1g", Nphoton_detect[3]);
        text.ReplaceAll("e+0","x10^{");
        text.Append("}");
        hCeren_d->SetTitle("Cerenkov photons ( N_{Reach_end}= "+ text+")");
        
    
        //normalize plots again
        hScinti->Scale(1.0 / NS_max);
        hScinti_d->Scale(1.0 / NS_det_max /2);
        hCeren->Scale(1.0 / NC_max);
        hCeren_d->Scale(1.0 / NC_det_max /2);
        hScinti->GetYaxis()->SetRangeUser(0, 2);
        hScinti->GetXaxis()->SetRangeUser(200, 1000);

        char title[100];
        if (is_bgo) sprintf(title, "1GeV muon (BGO, 2.4*2.4*10 cm^{3}, 2.4*2.4mm SiPM)");
        else sprintf(title, "1GeV muon (PWO, 2.4*2.4*10 cm^{3}, 2.4*2.4mm SiPM)");
        TPaveText *t = new TPaveText(0.35, 0.9, 0.6, 1.0, "brNDC");
        t->AddText(title);
        t->SetTextSize(0.05);
        t->SetBorderSize(0);
        t->SetFillColor(0);
        t->SetFillStyle(0);
        t->SetTextFont(42);
        t->Draw();
    
        //add self abosorption and PDEs
        double wc_low = 300;
        double wc_high = 1000;
        TF1 *RGB_sipm_QE = new TF1("RGB_sipm_QE", RGB_QE, 305.28, 900, 0);
        TF1 *UV_sipm_QE = new TF1("UV_sipm_QE", UV_QE, 361.161, 798.108, 0);
        RGB_sipm_QE->SetTitle("FBK RGB SiPM PDE");
        UV_sipm_QE->SetTitle("FBK NUV-HD SiPM PDE");
        RGB_sipm_QE->SetLineStyle(9);
        UV_sipm_QE->SetLineStyle(9);
        RGB_sipm_QE->SetLineWidth(3);
        UV_sipm_QE->SetLineWidth(3);
        RGB_sipm_QE->SetLineColor(46);
        UV_sipm_QE->SetLineColor(30);
        RGB_sipm_QE->Draw("same");
        UV_sipm_QE->Draw("same");
        TF1 *fABS_BGO;
        if (is_bgo)
            fABS_BGO = new TF1("fABS_BGO", ABS_BGO, 320, 590, 0);
        else
            fABS_BGO = new TF1("fABS_BGO", ABS_PWO, 367, 590, 0);
        fABS_BGO->SetTitle("Absorption coefficient");
        fABS_BGO->SetLineStyle(9);
        fABS_BGO->SetLineWidth(3);
        fABS_BGO->SetLineColor(9);
        fABS_BGO->Draw("same");
        auto leg = new TLegend(0.1, 0.9, 0.7, 0.9);
        leg = p1->BuildLegend();
        leg->SetLineWidth(0);
        leg->SetTextSize(0.03);
        leg->SetY1NDC(0.01);
        Style_t tfont = hScinti->GetYaxis()->GetTitleFont();
        Float_t tsize = hScinti->GetYaxis()->GetTitleSize();
        Style_t lfont = hScinti->GetYaxis()->GetLabelFont();
        Float_t lsize = hScinti->GetYaxis()->GetLabelSize();
        Double_t xmin = p1->GetUxmin();
        Double_t xmax = p1->GetUxmax();
        Double_t dx = (xmax - xmin) / 0.83; // 10 percent margins left and right
        Double_t ymin = hScinti->GetMinimum();
        Double_t ymax = hScinti->GetMaximum();
        TGaxis *axis = new TGaxis(1000, 0, 1000, 2, 0, 2, 510, "+L");
        axis->SetTitle("Quantum Effeciency");
        axis->SetTitleFont(tfont);
        axis->SetTitleSize(tsize);
        axis->SetTitleColor(kBlue);
        axis->SetTitleOffset(1.13);
        axis->SetLabelFont(lfont);
        axis->SetLabelSize(lsize);
        axis->SetLabelColor(kBlue);
        axis->SetLineColor(kBlue);
        axis->Draw();
        //delete c;
    }

    delete c1;
    // return numbers
    output_[0] = Nphoton_detect[0];
    output_[1] = Nphoton_detect[1];
    output_[2] = Nphoton_detect[2];
    output_[3] = Nphoton_detect[3];
    output_[4] = Efficiency;
    return true;
}

void make_angular_dependence()
{
    
    const int npoint = 4;  // number of points in angular scan
    double x_angles[npoint]={75,105,135,165}; //degree
    double y_C_S[npoint]; // C/S
    double N_Cherenkov[npoint]; // number of final detected C photons

    double output_[10];
    for (int i = 0; i < npoint; i++)
    {
        cout<<"open file "<<Form("test_%i_24_1GeVmuon.root", int(x_angles[i]))<<" ..."<<endl;
        // output_ will be [S_in_front, S_in_back, C_in_front, C_in_back, efficiency]
        // false -> pwo, true -> bgo;
        // false -> no plots, true -> make plots
        Get_Nphoton(Form("test_%i_24_1GeVmuon.root", int(x_angles[i])), output_, false, true); 
        x_angles[i] = x_angles[i] - 90;
        // Assume using the back for S, front for C
        y_C_S[i] = output_[2]/output_[1]; // C_in_front/S_in_back 
        N_Cherenkov[i] = output_[2] * output_[4] ; // C_in_front * efficiency
    }

    TCanvas *cfinal = new TCanvas();
    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    p1->SetGrid();
    p1->SetLeftMargin(0.14);
    p1->SetRightMargin(0.14);
    p1->Draw();
    p1->cd();

    TGraph *gcurve1 = new TGraph(npoint, x_angles, y_C_S);
    gcurve1->SetMarkerStyle(8);
    gcurve1->Draw("APL");
    gcurve1->SetMarkerColor(kBlue);

    // C/S plot
    TMultiGraph *gm = new TMultiGraph();
    gm->Add(gcurve1);
    gm->GetXaxis()->SetTitle("degree");
    gm->GetYaxis()->SetTitle("C/S");
    gm->GetXaxis()->SetLimits(-30, 120);
    gm->GetXaxis()->SetRangeUser(-30, 120);
    gm->GetYaxis()->SetLimits(0, 1);
    gm->GetYaxis()->SetRangeUser(0, 0.4);
    gm->Draw("AP");
 
    // N_Cherenkov Plot
    double norm_fc = 200;
    for (int i = 0; i < npoint; i++)
    {
        N_Cherenkov[i] /= norm_fc;
    }

    TGraph *g_NCheren = new TGraph(npoint, x_angles, N_Cherenkov);
    g_NCheren->SetMarkerStyle(4);
    g_NCheren->Draw("APL");
    g_NCheren->SetMarkerColor(kBlue);
    g_NCheren->SetLineColor(kBlue);

    gcurve1->SetTitle("C/S 24x24mm^{2}");
    g_NCheren->SetTitle("N_{C_photons} 24x24mm^{2}");
    gm->Add(g_NCheren);
    gm->Draw("AP");
    gPad->Update();

    Style_t tfont = gm->GetHistogram()->GetYaxis()->GetTitleFont();
    Float_t tsize = gm->GetHistogram()->GetYaxis()->GetTitleSize();
    Style_t lfont = gm->GetHistogram()->GetYaxis()->GetLabelFont();
    Float_t lsize = gm->GetHistogram()->GetYaxis()->GetLabelSize();
    Double_t xmin = p1->GetUxmin();
    Double_t xmax = p1->GetUxmax();
    Double_t dx = (xmax - xmin) / 0.83; // 10 percent margins left and right
    Double_t ymin = gm->GetHistogram()->GetMinimum();
    Double_t ymax = gm->GetHistogram()->GetMaximum();
    assert(ymin == 0);
    Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom

    TLegend *leg = gPad->BuildLegend();
    leg->SetFillColor(gPad->GetFillColor());
    leg->SetTextFont(lfont);
    leg->SetTextSize(0.75 * lsize);
    leg->Draw();
    gPad->Update();

    TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax * norm_fc, 510, "+L");
    axis->SetTitle("Number of detected C photons");
    axis->SetTitleFont(tfont);
    axis->SetTitleSize(tsize);
    axis->SetTitleColor(kBlue);
    axis->SetTitleOffset(1.13);
    axis->SetLabelFont(lfont);
    axis->SetLabelSize(lsize);
    axis->SetLabelColor(kBlue);
    axis->SetLineColor(kBlue);
    axis->Draw();

    TPaveText *t = new TPaveText(0.35, 0.9, 0.6, 1.0, "brNDC"); // left-up
    t->AddText(" PWO, Length 100 mm");
    t->SetTextSize(0.05);
    t->SetBorderSize(0);
    t->SetFillColor(0);
    t->SetFillStyle(0);
    t->SetTextFont(42);
    t->Draw();
}

void ana()
{
    make_angular_dependence();
}
