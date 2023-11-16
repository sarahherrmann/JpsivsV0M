#ifdef __CLING__
#pragma cling optimize(0)
#endif
void canvasResult()
{
//=========Macro generated from canvas: canvasResult/canvasResult
//=========  (Thu Nov  2 14:55:08 2023) by ROOT version 6.28/04
   TCanvas *canvasResult = new TCanvas("canvasResult", "canvasResult",0,53,800,600);
   gStyle->SetOptStat(0);
   canvasResult->Range(-0.8708919,-1.408696,7.870892,9.841304);
   canvasResult->SetFillColor(0);
   canvasResult->SetBorderMode(0);
   canvasResult->SetBorderSize(2);
   canvasResult->SetGridx();
   canvasResult->SetGridy();
   canvasResult->SetTopMargin(0.07478261);
   canvasResult->SetBottomMargin(0.1252174);
   canvasResult->SetFrameBorderMode(0);
   canvasResult->SetFrameBorderMode(0);
   
   TH1F *plot_copy__1 = new TH1F("plot_copy__1","",36,0,7);
   plot_copy__1->SetMinimum(0);
   plot_copy__1->SetMaximum(9);
   plot_copy__1->SetDirectory(nullptr);
   plot_copy__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   plot_copy__1->SetLineColor(ci);
   plot_copy__1->GetXaxis()->SetTitle("#frac{#it{N_{V0C}}}{<#it{N_{V0C}}>}");
   plot_copy__1->GetXaxis()->SetLabelFont(42);
   plot_copy__1->GetXaxis()->SetTitleOffset(1.2);
   plot_copy__1->GetXaxis()->SetTitleFont(42);
   plot_copy__1->GetYaxis()->SetTitle("#frac{dN_{J/#psi}/d#it{N_{V0C}}}{<dN_{J/#psi}/d#it{N_{V0C}}>}");
   plot_copy__1->GetYaxis()->SetLabelFont(42);
   plot_copy__1->GetYaxis()->SetTitleFont(42);
   plot_copy__1->GetZaxis()->SetLabelFont(42);
   plot_copy__1->GetZaxis()->SetTitleOffset(1);
   plot_copy__1->GetZaxis()->SetTitleFont(42);
   plot_copy__1->Draw("");
   
   Double_t JpsivsV0M_ptBin_0_fx1001[10] = { 0.195468, 0.5500184, 0.8514459, 1.120445, 1.478352, 1.837365, 2.16918, 2.659454, 3.490553, 5.196855 };
   Double_t JpsivsV0M_ptBin_0_fy1001[10] = { 0.3066861, 0.6304808, 0.9448298, 1.187347, 1.436202, 1.649669, 1.906703, 2.28576, 2.90412, 3.729401 };
   Double_t JpsivsV0M_ptBin_0_fex1001[10] = { 0.1917192, 0.1757426, 0.1182269, 0.1533754, 0.2108912, 0.1405941, 0.1981099, 0.3131414, 0.6390642, 0.9151249 };
   Double_t JpsivsV0M_ptBin_0_fey1001[10] = { 0.02627753, 0.02157655, 0.04454647, 0.05114904, 0.05448426, 0.06223935, 0.09352068, 0.07252045, 0.1210156, 0.2143924 };
   TGraphErrors *gre = new TGraphErrors(10,JpsivsV0M_ptBin_0_fx1001,JpsivsV0M_ptBin_0_fy1001,JpsivsV0M_ptBin_0_fex1001,JpsivsV0M_ptBin_0_fey1001);
   gre->SetName("JpsivsV0M_ptBin_0");
   gre->SetTitle("dN_{J/#psi}");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(8);
   
   TH1F *Graph_JpsivsV0M_ptBin_01001 = new TH1F("Graph_JpsivsV0M_ptBin_01001","dN_{J/#psi}",100,0,6.722803);
   Graph_JpsivsV0M_ptBin_01001->SetMinimum(0.2523678);
   Graph_JpsivsV0M_ptBin_01001->SetMaximum(4.310132);
   Graph_JpsivsV0M_ptBin_01001->SetDirectory(nullptr);
   Graph_JpsivsV0M_ptBin_01001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_JpsivsV0M_ptBin_01001->SetLineColor(ci);
   Graph_JpsivsV0M_ptBin_01001->GetXaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_01001->GetXaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_01001->GetXaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_01001->GetYaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_01001->GetYaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_01001->GetZaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_01001->GetZaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_01001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_JpsivsV0M_ptBin_01001);
   
   gre->Draw("p ");
   
   Double_t JpsivsV0M_ptBin_1_fx1002[10] = { 0.195468, 0.5500184, 0.8514459, 1.120445, 1.478352, 1.837365, 2.16918, 2.659454, 3.490553, 5.196855 };
   Double_t JpsivsV0M_ptBin_1_fy1002[10] = { 0.1986405, 0.5326775, 0.8191005, 1.103069, 1.448806, 1.773469, 2.187689, 2.690217, 3.537624, 4.692821 };
   Double_t JpsivsV0M_ptBin_1_fex1002[10] = { 0.1917192, 0.1757426, 0.1182269, 0.1533754, 0.2108912, 0.1405941, 0.1981099, 0.3131414, 0.6390642, 0.9151249 };
   Double_t JpsivsV0M_ptBin_1_fey1002[10] = { 0.007032919, 0.00640738, 0.0106795, 0.00310856, 0.003340086, 0.003116337, 0.03454157, 0.01085037, 0.006233284, 0.03018224 };
   gre = new TGraphErrors(10,JpsivsV0M_ptBin_1_fx1002,JpsivsV0M_ptBin_1_fy1002,JpsivsV0M_ptBin_1_fex1002,JpsivsV0M_ptBin_1_fey1002);
   gre->SetName("JpsivsV0M_ptBin_1");
   gre->SetTitle("dN_{J/#psi}");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(8);
   
   TH1F *Graph_JpsivsV0M_ptBin_11002 = new TH1F("Graph_JpsivsV0M_ptBin_11002","dN_{J/#psi}",100,0,6.722803);
   Graph_JpsivsV0M_ptBin_11002->SetMinimum(0.1724468);
   Graph_JpsivsV0M_ptBin_11002->SetMaximum(5.176143);
   Graph_JpsivsV0M_ptBin_11002->SetDirectory(nullptr);
   Graph_JpsivsV0M_ptBin_11002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_JpsivsV0M_ptBin_11002->SetLineColor(ci);
   Graph_JpsivsV0M_ptBin_11002->GetXaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_11002->GetXaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_11002->GetXaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_11002->GetYaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_11002->GetYaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_11002->GetZaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_11002->GetZaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_11002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_JpsivsV0M_ptBin_11002);
   
   gre->Draw("p ");
   
   Double_t JpsivsV0M_ptBin_2_fx1003[10] = { 0.195468, 0.5500184, 0.8514459, 1.120445, 1.478352, 1.837365, 2.16918, 2.659454, 3.490553, 5.196855 };
   Double_t JpsivsV0M_ptBin_2_fy1003[10] = { 0.1006378, 0.3592108, 0.6819162, 0.9744885, 1.456381, 1.946957, 2.434687, 3.18946, 4.269063, 6.806307 };
   Double_t JpsivsV0M_ptBin_2_fex1003[10] = { 0.1917192, 0.1757426, 0.1182269, 0.1533754, 0.2108912, 0.1405941, 0.1981099, 0.3131414, 0.6390642, 0.9151249 };
   Double_t JpsivsV0M_ptBin_2_fey1003[10] = { 0.003128911, 0.1419833, 0.003715379, 0.01769803, 0.01895322, 0.03612967, 0.04048521, 0.04815954, 0.05121334, 0.1542524 };
   gre = new TGraphErrors(10,JpsivsV0M_ptBin_2_fx1003,JpsivsV0M_ptBin_2_fy1003,JpsivsV0M_ptBin_2_fex1003,JpsivsV0M_ptBin_2_fey1003);
   gre->SetName("JpsivsV0M_ptBin_2");
   gre->SetTitle("dN_{J/#psi}");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(8);
   
   TH1F *Graph_JpsivsV0M_ptBin_21003 = new TH1F("Graph_JpsivsV0M_ptBin_21003","dN_{J/#psi}",100,0,6.722803);
   Graph_JpsivsV0M_ptBin_21003->SetMinimum(0.08775799);
   Graph_JpsivsV0M_ptBin_21003->SetMaximum(7.646865);
   Graph_JpsivsV0M_ptBin_21003->SetDirectory(nullptr);
   Graph_JpsivsV0M_ptBin_21003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_JpsivsV0M_ptBin_21003->SetLineColor(ci);
   Graph_JpsivsV0M_ptBin_21003->GetXaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_21003->GetXaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_21003->GetXaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_21003->GetYaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_21003->GetYaxis()->SetTitleFont(42);
   Graph_JpsivsV0M_ptBin_21003->GetZaxis()->SetLabelFont(42);
   Graph_JpsivsV0M_ptBin_21003->GetZaxis()->SetTitleOffset(1);
   Graph_JpsivsV0M_ptBin_21003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_JpsivsV0M_ptBin_21003);
   
   gre->Draw("p ");
   
   TLegend *leg = new TLegend(0.1453634,0.6747826,0.4260652,0.8747826,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("JpsivsV0M_ptBin_0","0.0 #leq #it{p_{T}} #leq 1.0 GeV/#it{c}","lpe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(8);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("JpsivsV0M_ptBin_1","1.0 #leq #it{p_{T}} #leq 3.0 GeV/#it{c}","lpe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(8);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("JpsivsV0M_ptBin_2","#it{p_{T}} > 3.0 GeV/#it{c}","lpe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(8);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   canvasResult->Modified();
   canvasResult->cd();
   canvasResult->SetSelected(canvasResult);
}
