#include <iostream>
#include <iomanip>
#include <fstream>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TGaxis.h>

#include "compAcc.h"

using namespace std;

void LoadAccFiles(string pfile1, string pfile2, string period)
{
  ifstream acc1_p(Form("%s/acceptance_pion_%s.txt",pfile1.c_str(),period.c_str()));
  ifstream acc1_k(Form("%s/acceptance_kaon_%s.txt",pfile1.c_str(),period.c_str()));
  ifstream acc1_pr(Form("%s/acceptance_proton_%s.txt",pfile1.c_str(),period.c_str()));
  ifstream acc1_h(Form("%s/acceptance_hadron_%s.txt",pfile1.c_str(),period.c_str()));
  ifstream acc2_p(Form("%s/acceptance_pion_%s.txt",pfile2.c_str(),period.c_str()));
  ifstream acc2_k(Form("%s/acceptance_kaon_%s.txt",pfile2.c_str(),period.c_str()));
  ifstream acc2_pr(Form("%s/acceptance_proton_%s.txt",pfile2.c_str(),period.c_str()));
  ifstream acc2_h(Form("%s/acceptance_hadron_%s.txt",pfile2.c_str(),period.c_str()));
  if (acc1_p.fail()) 
  {
    cerr << "\nERROR opening \""<< Form("%s/acceptance_pion_%s.txt",pfile1.c_str(),period.c_str()) <<"\": " << strerror(errno) << endl;
    abort();
  }

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
          for(int m = 0; m<2; m++)
          {
            acc1_p >> fAcceptance1[i][j][k].tab[1][m][0][0] >> fAcceptance1[i][j][k].tab[1][m][1][0];
            acc1_p >> fAcceptance1[i][j][k].tab[0][m][0][0] >> fAcceptance1[i][j][k].tab[0][m][1][0];
            acc1_k >> fAcceptance1[i][j][k].tab[1][m][0][1] >> fAcceptance1[i][j][k].tab[1][m][1][1];
            acc1_k >> fAcceptance1[i][j][k].tab[0][m][0][1] >> fAcceptance1[i][j][k].tab[0][m][1][1];
            acc1_pr >> fAcceptance1[i][j][k].tab[1][m][0][2] >> fAcceptance1[i][j][k].tab[1][m][1][2];
            acc1_pr >> fAcceptance1[i][j][k].tab[0][m][0][2] >> fAcceptance1[i][j][k].tab[0][m][1][2];
            acc1_h >> fAcceptance1[i][j][k].tab[1][m][0][3] >> fAcceptance1[i][j][k].tab[1][m][1][3];
            acc1_h >> fAcceptance1[i][j][k].tab[0][m][0][3] >> fAcceptance1[i][j][k].tab[0][m][1][3];
            acc2_p >> fAcceptance2[i][j][k].tab[1][m][0][0] >> fAcceptance2[i][j][k].tab[1][m][1][0];
            acc2_p >> fAcceptance2[i][j][k].tab[0][m][0][0] >> fAcceptance2[i][j][k].tab[0][m][1][0];
            acc2_k >> fAcceptance2[i][j][k].tab[1][m][0][1] >> fAcceptance2[i][j][k].tab[1][m][1][1];
            acc2_k >> fAcceptance2[i][j][k].tab[0][m][0][1] >> fAcceptance2[i][j][k].tab[0][m][1][1];
            acc2_pr >> fAcceptance2[i][j][k].tab[1][m][0][2] >> fAcceptance2[i][j][k].tab[1][m][1][2];
            acc2_pr >> fAcceptance2[i][j][k].tab[0][m][0][2] >> fAcceptance2[i][j][k].tab[0][m][1][2];
            acc2_h >> fAcceptance2[i][j][k].tab[1][m][0][3] >> fAcceptance2[i][j][k].tab[1][m][1][3];
            acc2_h >> fAcceptance2[i][j][k].tab[0][m][0][3] >> fAcceptance2[i][j][k].tab[0][m][1][3];
          }
      }
    }
  }
  acc1_p.close();
  acc1_k.close();
  acc1_pr.close();
  acc1_h.close();
  acc2_p.close();
  acc2_k.close();
  acc2_pr.close();
  acc2_h.close();
}

void yweightedavg()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int l=0; l<4; l++)
        {
            for(int m=0; m<2; m++)
            {
                for(int i=0; i<6; i++)
                {
                    if(fAcceptance1[x][i][z].tab[c][m][0][l])
                    {
                        fAcceptance1_yavg[x][z].tab[c][m][0][l]+=fAcceptance1[x][i][z].tab[c][m][0][l]/fAcceptance1[x][i][z].tab[c][m][1][l];
                        fAcceptance1_yavg[x][z].tab[c][m][1][l]+=1/fAcceptance1[x][i][z].tab[c][m][1][l];
                    }
                    if(fAcceptance2[x][i][z].tab[c][m][0][l])
                    {
                        fAcceptance2_yavg[x][z].tab[c][m][0][l]+=fAcceptance2[x][i][z].tab[c][m][0][l]/fAcceptance2[x][i][z].tab[c][m][1][l];
                        fAcceptance2_yavg[x][z].tab[c][m][1][l]+=1/fAcceptance2[x][i][z].tab[c][m][1][l];
                    }
                }
                if(fAcceptance1_yavg[x][z].tab[c][m][0][l])
                {
                    fAcceptance1_yavg[x][z].tab[c][m][1][l]=1/fAcceptance1_yavg[x][z].tab[c][m][1][l];
                    fAcceptance1_yavg[x][z].tab[c][m][0][l]*=fAcceptance1_yavg[x][z].tab[c][m][1][l];
                }
                if(fAcceptance2_yavg[x][z].tab[c][m][0][l])
                {
                    fAcceptance2_yavg[x][z].tab[c][m][1][l]=1/fAcceptance2_yavg[x][z].tab[c][m][1][l];
                    fAcceptance2_yavg[x][z].tab[c][m][0][l]*=fAcceptance2_yavg[x][z].tab[c][m][1][l];
                }
            }
        }
      }
    }
  }
}

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 2 *** Received : " << argc-1 << endl;
    cout << "./compAcc acc1 acc2 period" << endl;

    return 1;
  }
  
  LoadAccFiles(argv[1],argv[2],argv[3]);

  TCanvas c10("acc_comparison_hadron","acc_comparison_hadron",3200,1600);
  TCanvas c11("acc_comparison_pion","acc_comparison_pion",3200,1600);
  TCanvas c12("acc_comparison_kaon","acc_comparison_kaon",3200,1600);
  TCanvas c13("acc_comparison_proton","acc_comparison_proton",3200,1600);
  TCanvas c20("acc_comparison_yavg_hadron","acc_comparison_yavg_hadron",3200,1600);
  TCanvas c21("acc_comparison_yavg_pion","acc_comparison_yavg_pion",3200,1600);
  TCanvas c22("acc_comparison_yavg_kaon","acc_comparison_yavg_kaon",3200,1600);
  TCanvas c23("acc_comparison_yavg_proton","acc_comparison_yavg_proton",3200,1600);
  TLine l1(0.1,0.9,0.9,0.9);
  TLine l2(0.1,1.1,0.9,1.1);
  TLine l3(0.1,0.95,0.9,0.95);
  TLine l4(0.1,1.05,0.9,1.05);
  l3.SetLineStyle(2); l4.SetLineStyle(2);
  c10.SetFillColor(0);
  c11.SetFillColor(0);
  c12.SetFillColor(0);
  c13.SetFillColor(0);
  c20.SetFillColor(0);
  c21.SetFillColor(0);
  c22.SetFillColor(0);
  c23.SetFillColor(0);
  c10.Divide(9,5,0,0);
  c11.Divide(9,5,0,0);
  c12.Divide(9,5,0,0);
  c13.Divide(9,5,0,0);
  c20.Divide(5,2,0,0);
  c21.Divide(5,2,0,0);
  c22.Divide(5,2,0,0);
  c23.Divide(5,2,0,0);

  ofstream ofs_ra("ratiodiff.txt", std::ofstream::out | std::ofstream::trunc);

  TGraphErrors* R_h[2][2][9][6];
  TGraphErrors* R_p[2][2][9][6];
  TGraphErrors* R_k[2][2][9][6];
  TGraphErrors* R_pr[2][2][9][6];
  TGraphErrors* R_y_h[2][2][9];
  TGraphErrors* R_y_p[2][2][9];
  TGraphErrors* R_y_k[2][2][9];
  TGraphErrors* R_y_pr[2][2][9];

  for(int c=0; c<2; c++)
  {
    for(int m=0; m<2; m++)
    {
        for(int i=0; i<9; i++)
        {

        vector<double> r_y_h;
        vector<double> r_y_p;
        vector<double> r_y_k;
        vector<double> r_y_pr;
        vector<double> r_y_h_err;
        vector<double> r_y_p_err;
        vector<double> r_y_k_err;
        vector<double> r_y_pr_err;
        vector<double> z_range_r_y_p;
        vector<double> z_range_r_y_k;
        vector<double> z_range_r_y_pr;
        vector<double> z_range_r_y_h;

        for(int l=0; l<12; l++)
        {
            z_range_r_y_h.push_back(z_range[l]);
            z_range_r_y_p.push_back(z_range[l]);
            z_range_r_y_k.push_back(z_range[l]);
            z_range_r_y_pr.push_back(z_range[l]);
        }

        for(int j=0; j<6; j++)
        {
            std::vector<double> r_h;
            std::vector<double> r_p;
            std::vector<double> r_k;
            std::vector<double> r_pr;
            std::vector<double> r_h_err;
            std::vector<double> r_p_err;
            std::vector<double> r_k_err;
            std::vector<double> r_pr_err;
            std::vector<double> z_range_r_h;
            std::vector<double> z_range_r_p;
            std::vector<double> z_range_r_k;
            std::vector<double> z_range_r_pr;

            for(int l=0; l<12; l++)
            {
            z_range_r_h.push_back(z_range[l]);
            z_range_r_p.push_back(z_range[l]);
            z_range_r_k.push_back(z_range[l]);
            z_range_r_pr.push_back(z_range[l]);
            }

            for(int k=0; k<12; k++)
            {
            r_p.push_back(fAcceptance2[i][j][k].tab[c][m][0][0] ? fAcceptance1[i][j][k].tab[c][m][0][0]/fAcceptance2[i][j][k].tab[c][m][0][0] : 0);
            // cout << "fAcceptance2[" << i << "][" << j << "][" << k << "].tab[" << c << "][" << m << "][0][0] = " << fAcceptance2[i][j][k].tab[c][m][0][0] << endl;
            r_k.push_back(fAcceptance2[i][j][k].tab[c][m][0][1] ? fAcceptance1[i][j][k].tab[c][m][0][1]/fAcceptance2[i][j][k].tab[c][m][0][1] : 0);
            r_pr.push_back(fAcceptance2[i][j][k].tab[c][m][0][2] ? fAcceptance1[i][j][k].tab[c][m][0][2]/fAcceptance2[i][j][k].tab[c][m][0][2] : 0);
            r_h.push_back(fAcceptance2[i][j][k].tab[c][m][0][3] ? fAcceptance1[i][j][k].tab[c][m][0][3]/fAcceptance2[i][j][k].tab[c][m][0][3] : 0);
            r_p_err.push_back(sqrt((fAcceptance1[i][j][k].tab[c][m][1][0]+pow(fAcceptance2[i][j][k].tab[c][m][1][0],2)*fAcceptance1[i][j][k].tab[c][m][0][0]
                                    /pow(fAcceptance2[i][j][k].tab[c][m][0][0],2))/pow(fAcceptance2[i][j][k].tab[c][m][0][0],2)));
            r_k_err.push_back(sqrt((fAcceptance1[i][j][k].tab[c][m][1][1]+pow(fAcceptance2[i][j][k].tab[c][m][1][1],2)*fAcceptance1[i][j][k].tab[c][m][0][1]
                                    /pow(fAcceptance2[i][j][k].tab[c][m][0][1],2))/pow(fAcceptance2[i][j][k].tab[c][m][0][1],2)));
            r_pr_err.push_back(sqrt((fAcceptance1[i][j][k].tab[c][m][1][2]+pow(fAcceptance2[i][j][k].tab[c][m][1][2],2)*fAcceptance1[i][j][k].tab[c][m][0][2]
                                    /pow(fAcceptance2[i][j][k].tab[c][m][0][2],2))/pow(fAcceptance2[i][j][k].tab[c][m][0][2],2)));
            r_h_err.push_back(sqrt((fAcceptance1[i][j][k].tab[c][m][1][3]+pow(fAcceptance2[i][j][k].tab[c][m][1][3],2)*fAcceptance1[i][j][k].tab[c][m][0][3]
                                    /pow(fAcceptance2[i][j][k].tab[c][m][0][3],2))/pow(fAcceptance2[i][j][k].tab[c][m][0][3],2)));
            }

            for(int k=12; k>0; k--)
            {
            if(!r_h[k-1]) {r_h.erase(r_h.begin()+k-1); r_h_err.erase(r_h_err.begin()+k-1); z_range_r_h.erase(z_range_r_h.begin()+k-1);}
            if(!r_p[k-1]) {r_p.erase(r_p.begin()+k-1); r_p_err.erase(r_p_err.begin()+k-1); z_range_r_p.erase(z_range_r_p.begin()+k-1);}
            if(!r_k[k-1]) {r_k.erase(r_k.begin()+k-1); r_k_err.erase(r_k_err.begin()+k-1); z_range_r_k.erase(z_range_r_k.begin()+k-1);}
            if(!r_pr[k-1]) {r_pr.erase(r_pr.begin()+k-1); r_pr_err.erase(r_pr_err.begin()+k-1); z_range_r_pr.erase(z_range_r_pr.begin()+k-1);}
            }

            bool r_h_empty = 0;
            bool r_p_empty = 0;
            bool r_k_empty = 0;
            bool r_pr_empty = 0;

            if(!(int(r_h.size()))) {/* cout << "r_h is empty !" << endl; */ r_h_empty = 1;}
            // else cout << "r_h is not empty !" << endl;
            if(!(int(r_p.size()))) r_p_empty = 1;
            if(!(int(r_k.size()))) r_k_empty = 1;
            if(!(int(r_pr.size()))) r_pr_empty = 1;

            R_h[c][m][i][j] = new TGraphErrors(int(r_h.size()),&(z_range_r_h[0]),&(r_h[0]),0,&(r_h_err[0]));
            R_p[c][m][i][j] = new TGraphErrors(int(r_p.size()),&(z_range_r_p[0]),&(r_p[0]),0,&(r_p_err[0]));
            R_k[c][m][i][j] = new TGraphErrors(int(r_k.size()),&(z_range_r_k[0]),&(r_k[0]),0,&(r_k_err[0]));
            R_pr[c][m][i][j] = new TGraphErrors(int(r_pr.size()),&(z_range_r_pr[0]),&(r_pr[0]),0,&(r_pr_err[0]));

            if(!c && !m)
            {
            R_h[c][m][i][j]->SetMarkerColor(fMarkerColor[4]);
            R_p[c][m][i][j]->SetMarkerColor(fMarkerColor[4]);
            R_k[c][m][i][j]->SetMarkerColor(fMarkerColor[4]);
            R_pr[c][m][i][j]->SetMarkerColor(fMarkerColor[4]);
            }
            else if(!c && m)
            {
            R_h[c][m][i][j]->SetMarkerColor(fMarkerColor[3]);
            R_p[c][m][i][j]->SetMarkerColor(fMarkerColor[3]);
            R_k[c][m][i][j]->SetMarkerColor(fMarkerColor[3]);
            R_pr[c][m][i][j]->SetMarkerColor(fMarkerColor[3]);
            }
            else if(c && !m)
            {
            R_h[c][m][i][j]->SetMarkerColor(fMarkerColor[2]);
            R_p[c][m][i][j]->SetMarkerColor(fMarkerColor[2]);
            R_k[c][m][i][j]->SetMarkerColor(fMarkerColor[2]);
            R_pr[c][m][i][j]->SetMarkerColor(fMarkerColor[2]);
            }
            else
            {
            R_h[c][m][i][j]->SetMarkerColor(fMarkerColor[0]);
            R_p[c][m][i][j]->SetMarkerColor(fMarkerColor[0]);
            R_k[c][m][i][j]->SetMarkerColor(fMarkerColor[0]);
            R_pr[c][m][i][j]->SetMarkerColor(fMarkerColor[0]);
            }

            R_h[c][m][i][j]->SetMarkerSize(1);
            R_p[c][m][i][j]->SetMarkerSize(1);
            R_k[c][m][i][j]->SetMarkerSize(1);
            R_pr[c][m][i][j]->SetMarkerSize(1);

            R_h[c][m][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
            R_p[c][m][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
            R_k[c][m][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
            R_pr[c][m][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);

            if(!r_h_empty)
            {
            c10.cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(R_h[c][m][i][j])
            {
                if(!c)
                {
                R_h[c][m][i][j]->Draw("SAMEPA");
                R_h[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_h[c][m][i][j]->SetMinimum(0.);
                R_h[c][m][i][j]->SetMaximum(2.);
                R_h[c][m][i][j]->GetXaxis()->SetLabelSize(0.06);
                R_h[c][m][i][j]->GetYaxis()->SetLabelSize(0.06);
                R_h[c][m][i][j]->SetTitle("");
                if(j==4) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==4)
                {
                    R_h[c][m][i][j]->GetXaxis()->SetTitle("#font[12]{z}");
                    R_h[c][m][i][j]->GetXaxis()->SetTitleSize(0.08);
                    R_h[c][m][i][j]->GetXaxis()->SetTitleOffset(.8);
                }
                R_h[c][m][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_h[c][m][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0 && j==3)
                {
                    R_h[c][m][i][j]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{h}} #font[12]{ratio}");
                    R_h[c][m][i][j]->GetYaxis()->SetTitleSize(0.08);
                }
                R_h[c][m][i][j]->Draw("SAMEP");
                R_h[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_h[c][m][i][j]->SetMinimum(0.);
                R_h[c][m][i][j]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c10.Range(0.1,0.,0.9,2.);
                }
                else
                {
                R_h[c][m][i][j]->Draw("SAMEP");
                R_h[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_h[c][m][i][j]->SetMinimum(0.);
                R_h[c][m][i][j]->SetMaximum(2.);
                }
            }
            c10.Update();
            }

            if(!r_p_empty)
            {
            c11.cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(R_p[c][m][i][j])
            {
                if(!c)
                {
                R_p[c][m][i][j]->Draw("SAMEPA");
                R_p[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_p[c][m][i][j]->SetMinimum(0.);
                R_p[c][m][i][j]->SetMaximum(2.);
                R_p[c][m][i][j]->GetXaxis()->SetLabelSize(0.06);
                R_p[c][m][i][j]->GetYaxis()->SetLabelSize(0.06);
                R_p[c][m][i][j]->SetTitle("");
                if(j==4) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==4)
                {
                    R_p[c][m][i][j]->GetXaxis()->SetTitle("#font[12]{z}");
                    R_p[c][m][i][j]->GetXaxis()->SetTitleSize(0.08);
                    R_p[c][m][i][j]->GetXaxis()->SetTitleOffset(.8);
                }
                R_p[c][m][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_p[c][m][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0 && j==3)
                {
                    R_p[c][m][i][j]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{#pi}} #font[12]{ratio}");
                    R_p[c][m][i][j]->GetYaxis()->SetTitleSize(0.08);
                }
                R_p[c][m][i][j]->Draw("SAMEP");
                R_p[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_p[c][m][i][j]->SetMinimum(0.);
                R_p[c][m][i][j]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c11.Range(0.1,0.,0.9,2.);
                }
                else
                {
                R_p[c][m][i][j]->Draw("SAMEP");
                R_p[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_p[c][m][i][j]->SetMinimum(0.);
                R_p[c][m][i][j]->SetMaximum(2.);
                }
            }
            c11.Update();
            }

            if(!r_k_empty)
            {
            c12.cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(R_k[c][m][i][j])
            {
                if(!c)
                {
                R_k[c][m][i][j]->Draw("SAMEPA");
                R_k[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_k[c][m][i][j]->SetMinimum(0.);
                R_k[c][m][i][j]->SetMaximum(2.);
                R_k[c][m][i][j]->GetXaxis()->SetLabelSize(0.06);
                R_k[c][m][i][j]->GetYaxis()->SetLabelSize(0.06);
                R_k[c][m][i][j]->SetTitle("");
                if(j==4) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==4)
                {
                    R_k[c][m][i][j]->GetXaxis()->SetTitle("#font[12]{z}");
                    R_k[c][m][i][j]->GetXaxis()->SetTitleSize(0.08);
                    R_k[c][m][i][j]->GetXaxis()->SetTitleOffset(.8);
                }
                R_k[c][m][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_k[c][m][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0 && j==3)
                {
                    R_k[c][m][i][j]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{#pi}} #font[12]{ratio}");
                    R_k[c][m][i][j]->GetYaxis()->SetTitleSize(0.08);
                }
                R_k[c][m][i][j]->Draw("SAMEP");
                R_k[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_k[c][m][i][j]->SetMinimum(0.);
                R_k[c][m][i][j]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c12.Range(0.1,0.,0.9,2.);
                }
                else
                {
                R_k[c][m][i][j]->Draw("SAMEP");
                R_k[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_k[c][m][i][j]->SetMinimum(0.);
                R_k[c][m][i][j]->SetMaximum(2.);
                }
            }
            c12.Update();
            }

            if(!r_pr_empty)
            {
            c13.cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(R_pr[c][m][i][j])
            {
                if(!c)
                {
                R_pr[c][m][i][j]->Draw("SAMEPA");
                R_pr[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_pr[c][m][i][j]->SetMinimum(0.);
                R_pr[c][m][i][j]->SetMaximum(2.);
                R_pr[c][m][i][j]->GetXaxis()->SetLabelSize(0.06);
                R_pr[c][m][i][j]->GetYaxis()->SetLabelSize(0.06);
                R_pr[c][m][i][j]->SetTitle("");
                if(j==4) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==4)
                {
                    R_pr[c][m][i][j]->GetXaxis()->SetTitle("#font[12]{z}");
                    R_pr[c][m][i][j]->GetXaxis()->SetTitleSize(0.08);
                    R_pr[c][m][i][j]->GetXaxis()->SetTitleOffset(.8);
                }
                R_pr[c][m][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_pr[c][m][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0 && j==3)
                {
                    R_pr[c][m][i][j]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{p}} #font[12]{ratio}");
                    R_pr[c][m][i][j]->GetYaxis()->SetTitleSize(0.08);
                }
                R_pr[c][m][i][j]->Draw("SAMEP");
                R_pr[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_pr[c][m][i][j]->SetMinimum(0.);
                R_pr[c][m][i][j]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c13.Range(0.1,0.,0.9,2.);
                }
                else
                {
                R_pr[c][m][i][j]->Draw("SAMEP");
                R_pr[c][m][i][j]->GetXaxis()->SetLimits(0.1,0.9);
                R_pr[c][m][i][j]->SetMinimum(0.);
                R_pr[c][m][i][j]->SetMaximum(2.);
                }
            }
            c13.Update();
            }
        }

        yweightedavg();

        for(int k=0; k<12; k++)
        {
            double accy, accye;
            accy = fAcceptance2_yavg[i][k].tab[c][m][0][0] ? fAcceptance1_yavg[i][k].tab[c][m][0][0]/fAcceptance2_yavg[i][k].tab[c][m][0][0] : 0;
            accye = fAcceptance2_yavg[i][k].tab[c][m][0][0] ? sqrt((fAcceptance1_yavg[i][k].tab[c][m][1][0]+fAcceptance2_yavg[i][k].tab[c][m][1][0]*pow(fAcceptance1_yavg[i][k].tab[c][m][0][0],2)
                                    /pow(fAcceptance2_yavg[i][k].tab[c][m][0][0],2))/pow(fAcceptance2_yavg[i][k].tab[c][m][0][0],2)) : 0;
            r_y_p.push_back(accy);
            r_y_p_err.push_back(accye);
            int ratFlag;
            if(accy>=1)
            {
            if(accy+accye<=1) ratFlag = 1;
            else ratFlag = 0;
            }
            else
            {
            if(accy+accye>=1) ratFlag = 1;
            else ratFlag = 0;
            }
            ofs_ra << c << " " << i << " " << k << " " << accy << " " << accye << " " << ratFlag << endl;
        }

        for(int k=0; k<12; k++)
        {
            double accy, accye;
            accy = fAcceptance2_yavg[i][k].tab[c][m][0][1] ? fAcceptance1_yavg[i][k].tab[c][m][0][1]/fAcceptance2_yavg[i][k].tab[c][m][0][1] : 0;
            accye = fAcceptance2_yavg[i][k].tab[c][m][0][1] ? sqrt((fAcceptance1_yavg[i][k].tab[c][m][1][1]+fAcceptance2_yavg[i][k].tab[c][m][1][1]*pow(fAcceptance1_yavg[i][k].tab[c][m][0][1],2)
                                    /pow(fAcceptance2_yavg[i][k].tab[c][m][0][1],2))/pow(fAcceptance2_yavg[i][k].tab[c][m][0][1],2)) : 0;
            r_y_k.push_back(accy);
            r_y_k_err.push_back(accye);
            int ratFlag;
            if(accy>=1)
            {
            if(accy+accye<=1) ratFlag = 1;
            else ratFlag = 0;
            }
            else
            {
            if(accy+accye>=1) ratFlag = 1;
            else ratFlag = 0;
            }
            ofs_ra << c << " " << i << " " << k << " " << accy << " " << accye << " " << ratFlag << endl;
        }

        for(int k=0; k<12; k++)
        {
            double accy, accye;
            accy = fAcceptance2_yavg[i][k].tab[c][m][0][2] ? fAcceptance1_yavg[i][k].tab[c][m][0][2]/fAcceptance2_yavg[i][k].tab[c][m][0][2] : 0;
            accye = fAcceptance2_yavg[i][k].tab[c][m][0][2] ? sqrt((fAcceptance1_yavg[i][k].tab[c][m][1][2]+fAcceptance2_yavg[i][k].tab[c][m][1][2]*pow(fAcceptance1_yavg[i][k].tab[c][m][0][2],2)
                                    /pow(fAcceptance2_yavg[i][k].tab[c][m][0][2],2))/pow(fAcceptance2_yavg[i][k].tab[c][m][0][2],2)) : 0;
            r_y_pr.push_back(accy);
            r_y_pr_err.push_back(accye);
            int ratFlag;
            if(accy>=1)
            {
            if(accy+accye<=1) ratFlag = 1;
            else ratFlag = 0;
            }
            else
            {
            if(accy+accye>=1) ratFlag = 1;
            else ratFlag = 0;
            }
            ofs_ra << c << " " << i << " " << k << " " << accy << " " << accye << " " << ratFlag << endl;
        }

        for(int k=0; k<12; k++)
        {
            double accy, accye;
            accy = fAcceptance2_yavg[i][k].tab[c][m][0][3] ? fAcceptance1_yavg[i][k].tab[c][m][0][3]/fAcceptance2_yavg[i][k].tab[c][m][0][3] : 0;
            accye = fAcceptance2_yavg[i][k].tab[c][m][0][3] ? sqrt((fAcceptance1_yavg[i][k].tab[c][m][1][3]+fAcceptance2_yavg[i][k].tab[c][m][1][3]*pow(fAcceptance1_yavg[i][k].tab[c][m][0][3],2)
                                    /pow(fAcceptance2_yavg[i][k].tab[c][m][0][3],2))/pow(fAcceptance2_yavg[i][k].tab[c][m][0][3],2)) : 0;
            r_y_h.push_back(accy);
            r_y_h_err.push_back(accye);
            int ratFlag;
            if(accy>=1)
            {
            if(accy+accye<=1) ratFlag = 1;
            else ratFlag = 0;
            }
            else
            {
            if(accy+accye>=1) ratFlag = 1;
            else ratFlag = 0;
            }
            ofs_ra << c << " " << i << " " << k << " " << accy << " " << accye << " " << ratFlag << endl;
        }

        for(int k=12; k>0; k--)
        {
            if(!r_y_h[k-1]) {r_y_h.erase(r_y_h.begin()+k-1); r_y_h_err.erase(r_y_h_err.begin()+k-1); z_range_r_y_h.erase(z_range_r_y_h.begin()+k-1);}
            if(!r_y_p[k-1]) {r_y_p.erase(r_y_p.begin()+k-1); r_y_p_err.erase(r_y_p_err.begin()+k-1); z_range_r_y_p.erase(z_range_r_y_p.begin()+k-1);}
            if(!r_y_k[k-1]) {r_y_k.erase(r_y_k.begin()+k-1); r_y_k_err.erase(r_y_k_err.begin()+k-1); z_range_r_y_k.erase(z_range_r_y_k.begin()+k-1);}
            if(!r_y_pr[k-1]) {r_y_pr.erase(r_y_pr.begin()+k-1); r_y_pr_err.erase(r_y_pr_err.begin()+k-1); z_range_r_y_pr.erase(z_range_r_y_pr.begin()+k-1);}
        }

        bool r_y_h_empty = 0;
        bool r_y_p_empty = 0;
        bool r_y_k_empty = 0;
        bool r_y_pr_empty = 0;

        if(!(int(r_y_h.size()))) r_y_h_empty = 1;
        if(!(int(r_y_p.size()))) r_y_p_empty = 1;
        if(!(int(r_y_k.size()))) r_y_k_empty = 1;
        if(!(int(r_y_pr.size()))) r_y_pr_empty = 1;

        R_y_h[c][m][i] = new TGraphErrors(int(r_y_h.size()),&(z_range_r_y_h[0]),&(r_y_h[0]),0,&(r_y_h_err[0]));
        R_y_p[c][m][i] = new TGraphErrors(int(r_y_p.size()),&(z_range_r_y_p[0]),&(r_y_p[0]),0,&(r_y_p_err[0]));
        R_y_k[c][m][i] = new TGraphErrors(int(r_y_k.size()),&(z_range_r_y_k[0]),&(r_y_k[0]),0,&(r_y_k_err[0]));
        R_y_pr[c][m][i] = new TGraphErrors(int(r_y_pr.size()),&(z_range_r_y_pr[0]),&(r_y_pr[0]),0,&(r_y_pr_err[0]));

        if(!c && !m)
        {
            R_y_h[c][m][i]->SetMarkerColor(fMarkerColor[4]);
            R_y_p[c][m][i]->SetMarkerColor(fMarkerColor[4]);
            R_y_k[c][m][i]->SetMarkerColor(fMarkerColor[4]);
            R_y_pr[c][m][i]->SetMarkerColor(fMarkerColor[4]);
        }
        else if(!c && m)
        {
            R_y_h[c][m][i]->SetMarkerColor(fMarkerColor[3]);
            R_y_p[c][m][i]->SetMarkerColor(fMarkerColor[3]);
            R_y_k[c][m][i]->SetMarkerColor(fMarkerColor[3]);
            R_y_pr[c][m][i]->SetMarkerColor(fMarkerColor[3]);
        }
        else if(c && !m)
        {
            R_y_h[c][m][i]->SetMarkerColor(fMarkerColor[2]);
            R_y_p[c][m][i]->SetMarkerColor(fMarkerColor[2]);
            R_y_k[c][m][i]->SetMarkerColor(fMarkerColor[2]);
            R_y_pr[c][m][i]->SetMarkerColor(fMarkerColor[2]);
        }
        else
        {
            R_y_h[c][m][i]->SetMarkerColor(fMarkerColor[0]);
            R_y_p[c][m][i]->SetMarkerColor(fMarkerColor[0]);
            R_y_k[c][m][i]->SetMarkerColor(fMarkerColor[0]);
            R_y_pr[c][m][i]->SetMarkerColor(fMarkerColor[0]);
        }

        R_y_h[c][m][i]->SetMarkerSize(3);
        R_y_p[c][m][i]->SetMarkerSize(3);
        R_y_k[c][m][i]->SetMarkerSize(3);
        R_y_pr[c][m][i]->SetMarkerSize(3);

        R_y_h[c][m][i]->SetMarkerStyle(fMarkerStyle[0][c]);
        R_y_p[c][m][i]->SetMarkerStyle(fMarkerStyle[0][c]);
        R_y_k[c][m][i]->SetMarkerStyle(fMarkerStyle[0][c]);
        R_y_pr[c][m][i]->SetMarkerStyle(fMarkerStyle[0][c]);

        if(!r_y_h_empty)
        {
            c20.cd(i+1);
            gPad->SetFillStyle(4000);
            if(R_y_h[c][m][i])
            {
            if(!c)
            {
                R_y_h[c][m][i]->Draw("SAMEPA");
                R_y_h[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_h[c][m][i]->SetMinimum(0.);
                R_y_h[c][m][i]->SetMaximum(2.);
                R_y_h[c][m][i]->GetXaxis()->SetLabelSize(0.06);
                R_y_h[c][m][i]->GetYaxis()->SetLabelSize(0.06);
                R_y_h[c][m][i]->SetTitle("");
                if(i>4) gPad->SetBottomMargin(.15);
                if(i==0 || i==5) gPad->SetLeftMargin(.22);
                if(i==8)
                {
                R_y_h[c][m][i]->GetXaxis()->SetTitle("#font[12]{z}");
                R_y_h[c][m][i]->GetXaxis()->SetTitleSize(0.08);
                R_y_h[c][m][i]->GetXaxis()->SetTitleOffset(.8);
                }
                R_y_h[c][m][i]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_y_h[c][m][i]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0)
                {
                R_y_h[c][m][i]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{h}}#font[12]{ ratio}");
                R_y_h[c][m][i]->GetYaxis()->SetTitleSize(0.08);
                }
                R_y_h[c][m][i]->Draw("SAMEP");
                R_y_h[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_h[c][m][i]->SetMinimum(0.);
                R_y_h[c][m][i]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c20.Range(0.1,0.,0.9,2.);
            }
            else
            {
                R_y_h[c][m][i]->Draw("SAMEP");
                R_y_h[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_h[c][m][i]->SetMinimum(0.);
                R_y_h[c][m][i]->SetMaximum(2.0);
            }
            if(i==0)
            {
                c20.cd(10);
                if (!c && !m)
                {
                TLegend leg(0.1,0.8,0.9,1.);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{h^{-}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(!c && m)
                {
                TLegend leg(0.1,0.6,0.9,0.8);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{h^{-}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(c && !m)
                {
                TLegend leg(0.1,0.4,0.9,0.6);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{h^{+}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else
                {
                TLegend leg(0.1,0.2,0.9,0.4);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{h^{+}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                c20.cd(i+1);
            }
            }
            c20.Update();
        }

        if(!r_y_p_empty)
        {
            c21.cd(i+1);
            gPad->SetFillStyle(4000);
            if(R_y_p[c][m][i])
            {
            if(!c)
            {
                R_y_p[c][m][i]->Draw("SAMEPA");
                R_y_p[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_p[c][m][i]->SetMinimum(0.);
                R_y_p[c][m][i]->SetMaximum(2.);
                R_y_p[c][m][i]->GetXaxis()->SetLabelSize(0.06);
                R_y_p[c][m][i]->GetYaxis()->SetLabelSize(0.06);
                R_y_p[c][m][i]->SetTitle("");
                if(i>4) gPad->SetBottomMargin(.15);
                if(i==0 || i==5) gPad->SetLeftMargin(.22);
                if(i==8)
                {
                R_y_p[c][m][i]->GetXaxis()->SetTitle("#font[12]{z}");
                R_y_p[c][m][i]->GetXaxis()->SetTitleSize(0.08);
                R_y_p[c][m][i]->GetXaxis()->SetTitleOffset(.8);
                }
                R_y_p[c][m][i]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_y_p[c][m][i]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0)
                {
                R_y_p[c][m][i]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{#pi}}#font[12]{ ratio}");
                R_y_p[c][m][i]->GetYaxis()->SetTitleSize(0.08);
                }
                R_y_p[c][m][i]->Draw("SAMEP");
                R_y_p[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_p[c][m][i]->SetMinimum(0.);
                R_y_p[c][m][i]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c21.Range(0.1,0.,0.9,2.);
            }
            else
            {
                R_y_p[c][m][i]->Draw("SAMEP");
                R_y_p[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_p[c][m][i]->SetMinimum(0.);
                R_y_p[c][m][i]->SetMaximum(2.0);
            }
            if(i==0)
            {
                c21.cd(10);
                if (!c && !m)
                {
                TLegend leg(0.1,0.8,0.9,1.);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#pi^{-}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(!c && m)
                {
                TLegend leg(0.1,0.6,0.9,0.8);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#pi^{-}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(c && !m)
                {
                TLegend leg(0.1,0.4,0.9,0.6);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#pi^{+}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else
                {
                TLegend leg(0.1,0.2,0.9,0.4);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#pi^{+}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                c21.cd(i+1);
            }
            }
            c21.Update();
        }

        if(!r_y_k_empty)
        {
            c22.cd(i+1);
            gPad->SetFillStyle(4000);
            if(R_y_k[c][m][i])
            {
            if(!c)
            {
                R_y_k[c][m][i]->Draw("SAMEPA");
                R_y_k[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_k[c][m][i]->SetMinimum(0.);
                R_y_k[c][m][i]->SetMaximum(2.);
                R_y_k[c][m][i]->GetXaxis()->SetLabelSize(0.06);
                R_y_k[c][m][i]->GetYaxis()->SetLabelSize(0.06);
                R_y_k[c][m][i]->SetTitle("");
                if(i>4) gPad->SetBottomMargin(.15);
                if(i==0 || i==5) gPad->SetLeftMargin(.22);
                if(i==8)
                {
                R_y_k[c][m][i]->GetXaxis()->SetTitle("#font[12]{z}");
                R_y_k[c][m][i]->GetXaxis()->SetTitleSize(0.08);
                R_y_k[c][m][i]->GetXaxis()->SetTitleOffset(.8);
                }
                R_y_k[c][m][i]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_y_k[c][m][i]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0)
                {
                R_y_k[c][m][i]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{K}}#font[12]{ ratio}");
                R_y_k[c][m][i]->GetYaxis()->SetTitleSize(0.08);
                }
                R_y_k[c][m][i]->Draw("SAMEP");
                R_y_k[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_k[c][m][i]->SetMinimum(0.);
                R_y_k[c][m][i]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c22.Range(0.1,0.,0.9,2.);
            }
            else
            {
                R_y_k[c][m][i]->Draw("SAMEP");
                R_y_k[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_k[c][m][i]->SetMinimum(0.);
                R_y_k[c][m][i]->SetMaximum(2.0);
            }
            if(i==0)
            {
                c22.cd(10);
                if (!c && !m)
                {
                TLegend leg(0.1,0.8,0.9,1.);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{K^{-}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(!c && m)
                {
                TLegend leg(0.1,0.6,0.9,0.8);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{K^{-}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(c && !m)
                {
                TLegend leg(0.1,0.4,0.9,0.6);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{K^{+}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else
                {
                TLegend leg(0.1,0.2,0.9,0.4);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{K^{+}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                c22.cd(i+1);
            }
            }
            c22.Update();
        }

        if(!r_y_pr_empty)
        {
            c23.cd(i+1);
            gPad->SetFillStyle(4000);
            if(R_y_pr[c][m][i])
            {
            if(!c)
            {
                R_y_pr[c][m][i]->Draw("SAMEPA");
                R_y_pr[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_pr[c][m][i]->SetMinimum(0.);
                R_y_pr[c][m][i]->SetMaximum(2.);
                R_y_pr[c][m][i]->GetXaxis()->SetLabelSize(0.06);
                R_y_pr[c][m][i]->GetYaxis()->SetLabelSize(0.06);
                R_y_pr[c][m][i]->SetTitle("");
                if(i>4) gPad->SetBottomMargin(.15);
                if(i==0 || i==5) gPad->SetLeftMargin(.22);
                if(i==8)
                {
                R_y_pr[c][m][i]->GetXaxis()->SetTitle("#font[12]{z}");
                R_y_pr[c][m][i]->GetXaxis()->SetTitleSize(0.08);
                R_y_pr[c][m][i]->GetXaxis()->SetTitleOffset(.8);
                }
                R_y_pr[c][m][i]->GetXaxis()->SetNdivisions(304,kTRUE);
                R_y_pr[c][m][i]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==0)
                {
                R_y_pr[c][m][i]->GetYaxis()->SetTitle("#font[12]{Acc}^{#font[12]{p}}#font[12]{ ratio}");
                R_y_pr[c][m][i]->GetYaxis()->SetTitleSize(0.08);
                }
                R_y_pr[c][m][i]->Draw("SAMEP");
                R_y_pr[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_pr[c][m][i]->SetMinimum(0.);
                R_y_pr[c][m][i]->SetMaximum(2.);
                l1.Draw("SAME");
                l2.Draw("SAME");
                l3.Draw("SAME");
                l4.Draw("SAME");
                c23.Range(0.1,0.,0.9,2.);
            }
            else
            {
                R_y_pr[c][m][i]->Draw("SAMEP");
                R_y_pr[c][m][i]->GetXaxis()->SetLimits(0.1,0.9);
                R_y_pr[c][m][i]->SetMinimum(0.);
                R_y_pr[c][m][i]->SetMaximum(2.0);
            }
            if(i==0)
            {
                c23.cd(10);
                if (!c && !m)
                {
                TLegend leg(0.1,0.8,0.9,1.);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#bar{p}; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(!c && m)
                {
                TLegend leg(0.1,0.6,0.9,0.8);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{#bar{p}; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else if(c && !m)
                {
                TLegend leg(0.1,0.4,0.9,0.6);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{p; #mu^{-} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                else
                {
                TLegend leg(0.1,0.2,0.9,0.4);
                leg.SetFillColor(0);
                leg.AddEntry(R_y_h[c][m][i],"#font[12]{p; #mu^{+} beam}");
                leg.SetBorderSize(0);
                leg.DrawClone("Same");
                }
                c23.cd(i+1);
            }
            }
            c23.Update();
        }
        }
    }
  }

  TLatex fTitle;

  c20.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c20.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c20.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c20.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c20.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c20.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c20.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c20.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c20.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c21.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c21.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c21.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c21.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c21.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c21.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c21.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c21.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c21.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c22.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c22.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c22.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c22.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c22.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c22.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c22.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c22.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c22.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c23.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c23.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c23.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c23.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c23.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c23.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c23.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c23.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c23.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c20.Update();
  c21.Update();
  c22.Update();
  c23.Update();

  c10.Print("compAcc/acc_ratio.pdf(");
  c11.Print("compAcc/acc_ratio.pdf");
  c12.Print("compAcc/acc_ratio.pdf");
  c13.Print("compAcc/acc_ratio.pdf)");
  c20.Print("compAcc/acc_ratio_yavg.pdf(");
  c21.Print("compAcc/acc_ratio_yavg.pdf");
  c22.Print("compAcc/acc_ratio_yavg.pdf");
  c23.Print("compAcc/acc_ratio_yavg.pdf)");

  ofs_ra.close();

  return 0;
}
