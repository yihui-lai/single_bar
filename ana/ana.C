#include "TStyle.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include <TMath.h>
#include <String.h>
#include <iostream>
#include <utility>

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

void C_S(TH1F *hc, TH1F *hs)
{
    //plot C ratio vs S ratio from 300 nm to 1000 nm
    TCanvas *cratio = new TCanvas();
    hc->Scale(1.0 / hc->Integral(300, 1000));
    hs->Scale(1.0 / hs->Integral(300, 1000));
    const int np = 351;
    double cc[np], ss[np];
    for (int i = 0; i < np; i++)
    {
        cc[i] = hc->Integral(300, 300 + 2 * i);
        ss[i] = hs->Integral(300, 300 + 2 * i);
        //if(hs->Integral(300,300+2*i)!=0 && (ss[i]-ss[i-1])!=0 && (cc[i]-cc[i-1])/(ss[i]-ss[i-1])>0.99 && i>1) cout<<300+2*i<<" "<<(cc[i]-cc[i-1])/(ss[i]-ss[i-1])<<endl;
        if (300 + 2 * i == 650)
            cout << 300 + 2 * i << " " << (ss[i]) << " " << (cc[i]) << endl;
        if (300 + 2 * i == 380)
            cout << 300 + 2 * i << " " << (ss[i]) << " " << (cc[i]) << endl;
        //else cc[i] =0;
        //ss[i] = 300+7.0*i;
    }
    TGraph *cscurve = new TGraph(np, ss, cc);
    cscurve->SetLineColor(kBlue);
    cscurve->SetLineWidth(2);
    cscurve->SetTitle("; N_{Scinti}/N_{Scinti_tot}; N_{Ceren}/N_{Ceren_tot}");
    cscurve->Draw("APL");
    TPaveText *t = new TPaveText(0.35, 0.9, 0.6, 1.0, "brNDC");
    t->SetTextSize(0.05);
    //t->AddText("PbWO_{4}, 1*1*10 cm^{3}");
    //t->AddText("BGO, 4*4*4 cm^{3}");
    t->SetBorderSize(0);
    t->SetFillColor(0);
    t->SetFillStyle(0);
    t->SetTextFont(42);
    t->Draw();
}

double Getratio(TH1F *hdetection, bool use_UV, double wvmin, double wvmax)
{
    //get the effective detection effeciency
    assert(wvmin <= wvmax);
    //hdetection->Draw();
    double wc_low = wvmin;
    double wc_high = wvmax;
    TF1 *rat_RGB_sipm_QE = new TF1("rat_RGB_sipm_QE", RGB_QE, wc_low, wc_high, 0);
    TF1 *rat_UV_sipm_QE = new TF1("rat_UV_sipm_QE", UV_QE, wc_low, wc_high, 0);
    double rat = 0;
    for (int ix = wc_low; ix < wc_high; ix++)
    {
        //cout<<"ix: "<<ix<<" Bincenter: "<< hdetection->GetBinCenter(ix) <<" Bincontent "<<hdetection->GetBinContent(ix)<<" sipm: "<<rat_RGB_sipm_QE->Eval(hdetection->GetBinCenter(ix))<<endl;
        if (use_UV)
        {
            rat += rat_UV_sipm_QE->Eval(hdetection->GetBinCenter(ix)) * hdetection->GetBinContent(ix);
        }
        else
        {
            rat += rat_RGB_sipm_QE->Eval(hdetection->GetBinCenter(ix)) * hdetection->GetBinContent(ix);
        }
    }
    //cout<<"ratio: "<<rat/hdetection->Integral(300,1000)<<endl;
    return rat / hdetection->Integral(300, 1000);
}

double Get_Nphoton(char *filename, double *N_C, bool is_bgo)
{
    //Get the number of photons
    double hh[6];
    int hhh[6];
    TFile *f = new TFile(filename);
    char fname[100];
    if (is_bgo)
        sprintf(fname, "1GeV muon (BGO, 1*1*10 cm^{3}, 6*6mm SiPM)");
    else
        sprintf(fname, "1GeV muon (PWO, 1*1*10 cm^{3}, 6*6mm SiPM)");

    TCanvas *c1 = new TCanvas();
    TH1F *hcount5 = new TH1F("hcount5", "", 1000, 0, 1000000);
    ((TTree *)f->Get("tree"))->Draw("ECAL_f_total_S>>hcount5");
    //TCanvas* c2=new TCanvas();
    TH1F *hcount6 = new TH1F("hcount6", "", 1000, 0, 200000);
    ((TTree *)f->Get("tree"))->Draw("ECAL_f_total_C>>hcount6");
    hh[4] = hcount5->GetMean();
    hh[5] = hcount6->GetMean();
    cout << "tot S:" << hh[4] << " tot C: " << hh[5] << endl;

    //TFile* f=new TFile("elec_9MeV_BGO.root");
    //TCanvas* c3=new TCanvas();
    TH1F *hcount1 = new TH1F("hcount1", "", 5000, 0, hh[4] / 2);
    ((TTree *)f->Get("tree"))->Draw("SDdetected_ff_S>>hcount1");

    cout << "S hit front " << hcount1->GetMean() << endl;
    TH1F *hcount2 = new TH1F("hcount2", "", 5000, 0, hh[4] / 2);
    //TCanvas* c4=new TCanvas();
    ((TTree *)f->Get("tree"))->Draw("SDdetected_rr_S>>hcount2");
    cout << "S hit back " << hcount2->GetMean() << endl;
    TH1F *hcount3 = new TH1F("hcount3", "", 5000, 0, hh[5] / 2);
    //TCanvas* c5=new TCanvas();
    ((TTree *)f->Get("tree"))->Draw("SDdetected_ff_C>>hcount3");
    cout << "C hit front " << hcount3->GetMean() << endl;
    TH1F *hcount4 = new TH1F("hcount4", "", 5000, 0, hh[5] / 2);
    //TCanvas* c6=new TCanvas();
    ((TTree *)f->Get("tree"))->Draw("SDdetected_rr_C>>hcount4");
    cout << "C hit back " << hcount4->GetMean() << endl;
    //TCanvas* c7=new TCanvas();
    hh[0] = hcount1->GetMean();
    hh[1] = hcount2->GetMean();
    hh[2] = hcount3->GetMean();
    hh[3] = hcount4->GetMean();

    for (int i = 0; i < 6; i++)
    {
        hhh[i] = 0;
        if (hh[i] >= 1 && hh[i] < 10)
        {
            hh[i] = hh[i];
            hhh[i] = 0;
        }
        else if (hh[i] >= 10 && hh[i] < 100)
        {
            hh[i] = hh[i] / 10;
            hhh[i] = 1;
        }
        else if (hh[i] >= 100 && hh[i] < 1000)
        {
            hh[i] = hh[i] / 100;
            hhh[i] = 2;
        }
        else if (hh[i] < 10000)
        {
            hh[i] = hh[i] / 1000;
            hhh[i] = 3;
        }
        else if (hh[i] < 100000)
        {
            hh[i] = hh[i] / 10000;
            hhh[i] = 4;
        }
        else if (hh[i] < 1000000)
        {
            hh[i] = hh[i] / 100000;
            hhh[i] = 5;
        }
        else if (hh[i] < 10000000)
        {
            hh[i] = hh[i] / 1000000;
            hhh[i] = 6;
        }
        else if (hh[i] < 100000000)
        {
            hh[i] = hh[i] / 10000000;
            hhh[i] = 7;
        }
        else
        {
            cout << "too large" << endl;
            return 0;
        }
    }

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

    //TFile* f=new TFile("muon_1GeV_silicone.root");
    //TFile* f=new TFile("elec_9MeV_silicone.root");
    int num = ((TTree *)f->Get("tree"))->GetEntries();
    TH1F *hCeren = (TH1F *)f->Get("h_phot_lambda_ECAL_f_produce_Ceren");
    hCeren->Scale(1.0 / num);
    TH1F *hScinti = (TH1F *)f->Get("h_phot_lambda_ECAL_f_produce_Scin");
    hScinti->Scale(1.0 / num);

    double sf1 = hScinti->GetMaximum();
    double sf2 = hCeren->GetMaximum();
    hCeren->SetLineColor(kRed);
    hScinti->SetLineColor(kBlue);
    hCeren->SetLineWidth(3);
    hScinti->SetLineWidth(3);
    hScinti->GetXaxis()->SetTitle("Wavelength (nm)");
    hScinti->GetYaxis()->SetTitle("Number of photons (normalized)");
    hScinti->Draw("HIST");
    hCeren->Draw("same HIST");
    //cout<<"N_gen: S,"<<hScinti->Integral()<<" C,"<<hCeren->Integral()<<endl;

    TH1F *hCeren_d = (TH1F *)f->Get("h_phot_lambda_ECAL_f_Ceren");
    hCeren_d->Scale(1.0 / num);
    TH1F *hScinti_d = (TH1F *)f->Get("h_phot_lambda_ECAL_f_Scin");
    hScinti_d->Scale(1.0 / num);
    double sfd1 = hScinti_d->GetMaximum();
    double sfd2 = hCeren_d->GetMaximum();

    hCeren_d->SetLineColor(kBlack);
    hScinti_d->SetLineColor(kOrange);
    hCeren_d->SetLineWidth(3);
    hScinti_d->SetLineWidth(3);
    hScinti_d->Draw("same HIST");
    hCeren_d->Draw("same HIST");
    //cout<<"N_detect: S,"<<hScinti_d->Integral()<<" C,"<<hCeren_d->Integral()<<endl;

    hScinti->SetTitle(Form("Scintillation photons ( N_{Generated}= %.1f*10^{%i} )", hh[4], hhh[4]));
    hCeren->SetTitle(Form("Cerenkov photons ( N_{Generated}= %.1f*10^{%i} )", hh[5], hhh[5]));
    hScinti_d->SetTitle(Form("Scintillation photons at one end  ( N_{Reach_end}= %.1f*10^{%i} )", hh[0], hhh[0]));
    hCeren_d->SetTitle(Form("Cerenkov photons at one end  ( N_{Reach_end}= %.1f*10^{%i} )", hh[3], hhh[3]));

    //cout<<"max:"<<sf1<<" , "<<sf2<<endl;

    hScinti->Scale(1.0 / sf1);
    hScinti_d->Scale(1.0 / sfd1 / 2);
    hCeren->Scale(1.0 / sf2);
    hCeren_d->Scale(1.0 / sfd2 / 2);
    hScinti->GetYaxis()->SetRangeUser(0, 2);
    hScinti->GetXaxis()->SetRangeUser(200, 1000);

    TPaveText *t = new TPaveText(0.35, 0.9, 0.6, 1.0, "brNDC"); // left-up
    //t->AddText(" 9MeV electron (BGO, 4*4*4 cm^{3})");
    t->AddText(fname);
    t->SetTextSize(0.05);
    //t->AddText("1GeV muon (PbWO_{4}, 1*1*10 cm^{3})");
    //t->AddText(" 9MeV electron (PbWO_{4}, 1*1*10 cm^{3})");
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

    //TF1 *fABS_BGO = new TF1("fABS_BGO", ABS_BGO, 320, 590, 0);
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
    //axis->SetTitle("Absorption coefficient (cm^{-1})");
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

    N_C[0] = hh[0] * pow(10, hhh[0]);
    N_C[1] = hh[1] * pow(10, hhh[1]);
    N_C[2] = hh[2] * pow(10, hhh[2]);
    N_C[3] = hh[3] * pow(10, hhh[3]);

    N_C[4] = Getratio((TH1F *)f->Get("h_phot_lambda_ECAL_f_Ceren"), false, 300, 380) + Getratio((TH1F *)f->Get("h_phot_lambda_ECAL_f_Ceren"), false, 650, 1000);

    delete c1;
    //delete c;
    return N_C[2] / N_C[1];
}

void ana_angle_dependence()
{

    const int np = 9;
    double xx[np];
    double curve3[np], curve4[np];
    double curve1[np], curve2[np];

    double NC40[np], NC120[np], NC200[np], NC280[np];

    double N_C[10];
    /*
    for(int i=0;i<np;i++){
        xx[i]=i*10+75-90;
        Get_Nphoton(Form("bgo_20_20_200//muon_1GeV_N20_BGO_%idegree_40mm.root",i*10+75), N_C,true);
        curve1[i]=N_C[2]/N_C[1];
        NC40[i]=N_C[2]*N_C[4];
        Get_Nphoton(Form("bgo_20_20_200//muon_1GeV_N20_BGO_%idegree_120mm.root",i*10+75), N_C,true);
        curve2[i]=N_C[2]/N_C[1];
        NC120[i]=N_C[2]*N_C[4];
        Get_Nphoton(Form("bgo_20_20_200//muon_1GeV_N20_BGO_%idegree_200mm.root",i*10+75), N_C,true);
        curve3[i]=N_C[2]/N_C[1];
        NC200[i]=N_C[2]*N_C[4];
        Get_Nphoton(Form("bgo_20_20_200//muon_1GeV_N20_BGO_%idegree_280mm.root",i*10+75), N_C,true);
        curve4[i]=N_C[2]/N_C[1];
        NC280[i]=N_C[2]*N_C[4];
    }
    */
    for (int i = 0; i < np; i++)
    {
        xx[i] = i * 10 + 75 - 90;
        curve1[i] = Get_Nphoton(Form("angle_plot/muon_1GeV_N50_BGO_%idegree_500mm_A10mm_sipm6mm.root", i * 10 + 75), N_C, true);
        NC40[i] = N_C[2] * N_C[4];
        cout << "detection coeff: " << N_C[4] << endl;
        curve2[i] = Get_Nphoton(Form("angle_plot/muon_1GeV_N50_BGO_%idegree_500mm_A50mm_sipm6mm.root", i * 10 + 75), N_C, true);
        NC120[i] = N_C[2] * N_C[4];
        cout << "detection coeff: " << N_C[4] << endl;
        curve3[i] = Get_Nphoton(Form("angle_plot/muon_1GeV_N50_BGO_%idegree_500mm_A100mm_sipm6mm.root", i * 10 + 75), N_C, true);
        NC200[i] = N_C[2] * N_C[4];
        cout << "detection coeff: " << N_C[4] << endl;
        curve4[i] = Get_Nphoton(Form("angle_plot/muon_1GeV_N50_BGO_%idegree_500mm_A150mm_sipm6mm.root", i * 10 + 75), N_C, true);
        NC280[i] = N_C[2] * N_C[4];
        cout << "detection coeff: " << N_C[4] << endl;
    }

    TGraph *gcurve1 = new TGraph(np, xx, curve1);
    TGraph *gcurve2 = new TGraph(np, xx, curve2);
    TGraph *gcurve3 = new TGraph(np, xx, curve3);
    TGraph *gcurve4 = new TGraph(np, xx, curve4);
    gcurve1->SetMarkerStyle(8);
    gcurve1->Draw("APL");
    gcurve1->SetMarkerColor(kBlue);
    gcurve1->SetLineColor(kBlue);
    gcurve2->Draw("APL");
    gcurve2->SetMarkerStyle(8);
    gcurve2->SetMarkerColor(kRed);
    gcurve2->SetLineColor(kRed);

    gcurve3->SetMarkerStyle(22);
    gcurve3->Draw("APL");
    gcurve3->SetMarkerColor(kBlue);
    gcurve3->SetLineColor(kBlue);
    gcurve4->Draw("APL");
    gcurve4->SetMarkerStyle(22);
    gcurve4->SetMarkerColor(kRed);
    gcurve4->SetLineColor(kRed);

    TCanvas *c = new TCanvas();

    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    p1->SetGrid();
    p1->SetLeftMargin(0.14);
    p1->SetRightMargin(0.14);
    p1->Draw();
    p1->cd();

    TMultiGraph *gm = new TMultiGraph();
    gm->Add(gcurve1);
    gm->Add(gcurve2);
    gm->Add(gcurve3);
    gm->Add(gcurve4);
    gm->GetXaxis()->SetTitle("degree");
    gm->GetYaxis()->SetTitle("C/S");
    gm->GetXaxis()->SetLimits(-30, 120);
    gm->GetXaxis()->SetRangeUser(-30, 120);
    gm->GetYaxis()->SetLimits(0, 1);
    gm->GetYaxis()->SetRangeUser(0, 0.02);
    gm->Draw("AP");

    double norm_fc = 8000;
    for (int i = 0; i < np; i++)
    {
        NC40[i] /= norm_fc;
        NC120[i] /= norm_fc;
        NC200[i] /= norm_fc;
        NC280[i] /= norm_fc;
    }

    TGraph *gc40 = new TGraph(np, xx, NC40);
    TGraph *gc120 = new TGraph(np, xx, NC120);
    TGraph *gc200 = new TGraph(np, xx, NC200);
    TGraph *gc280 = new TGraph(np, xx, NC280);
    gc40->SetMarkerStyle(4);
    gc40->Draw("APL");
    gc40->SetMarkerColor(kBlue);
    gc40->SetLineColor(kBlue);
    gc120->Draw("APL");
    gc120->SetMarkerStyle(4);
    gc120->SetMarkerColor(kRed);
    gc120->SetLineColor(kRed);

    gc200->SetMarkerStyle(32);
    gc200->Draw("APL");
    gc200->SetMarkerColor(kBlue);
    gc200->SetLineColor(kBlue);
    gc280->Draw("APL");
    gc280->SetMarkerStyle(32);
    gc280->SetMarkerColor(kRed);
    gc280->SetLineColor(kRed);

    gcurve1->SetTitle("C/S 10x10mm^{2}");
    gcurve2->SetTitle("C/S 50x50mm^{2}");
    gcurve3->SetTitle("C/S 100x100mm^{2}");
    gcurve4->SetTitle("C/S 150x150mm^{2}");
    gc200->SetTitle("N_{C_photons} 100x100mm^{2}");
    gc280->SetTitle("N_{C_photons} 150x150mm^{2}");
    gc40->SetTitle("N_{C_photons} 10x10mm^{2}");
    gc120->SetTitle("N_{C_photons} 50x50mm^{2}");
    /*
     gcurve1->SetTitle("C/S 40mm");
     gcurve2->SetTitle("C/S 120mm");
     gcurve3->SetTitle("C/S 200mm");
     gcurve4->SetTitle("C/S 280mm");
    gc40->SetTitle("N_{C_photons} 40mm");
    gc120->SetTitle("N_{C_photons} 120mm");
    gc200->SetTitle("N_{C_photons} 200mm");
    gc280->SetTitle("N_{C_photons} 280mm");
*/
    gm->Add(gc40);
    gm->Add(gc120);
    gm->Add(gc200);
    gm->Add(gc280);
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

    //TLegend *leg = new TLegend(0.2, 0.77, 0.7, 0.90);
    leg->SetFillColor(gPad->GetFillColor());
    leg->SetTextFont(lfont);
    leg->SetTextSize(0.75 * lsize);
    leg->Draw();
    gPad->Update();

    //TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
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
    //t->AddText(" BGO, area 20*20 mm^{2}");
    t->AddText(" BGO, Length 500 mm");
    t->SetTextSize(0.05);
    t->SetBorderSize(0);
    t->SetFillColor(0);
    t->SetFillStyle(0);
    t->SetTextFont(42);
    t->Draw();
}

void ana()
{
    //ana_angle_dependence();

    double N_C[10];
    Get_Nphoton("bgo_10_10_100//muon_1GeV_N20_BGO_90degree_100mm_A10mm_sipm6mm.root", N_C, true);
    //Get_Nphoton("realistic/muon_1GeV_N50_BGO_90degree_40mm_A40mm_sipm6mm.root", N_C, true);
    //TFile* ff=new TFile("bgo_10_10_100//muon_1GeV_N20_BGO_90degree_100mm_A10mm_sipm6mm.root");
    //C_S((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Ceren"),(TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"));
    //cout<<"S in 300 to 380: "<<Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 300, 380)<<endl;
    //cout<<"S in 650 to 1000: "<<Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 650, 1000)<<endl;
    //cout<<"S in 380 to 650: "<<Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 380, 650)<<endl;
    
    /*
    Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 300, 380);
    Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 650, 1000);
    Getratio((TH1F*) ff->Get("h_phot_lambda_ECAL_f_Scin"), false, 380, 650);
    */
}
