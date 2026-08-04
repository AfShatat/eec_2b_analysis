#include "CMSStyle.C"

// usage:  root -l rootlogon.C plotNice_UParT_roc.C


void plotNice_UParT_roc(){

  setCMSStyle();

  auto *fin = new TFile("btag_roc_qcd.root");

  auto *g1b = (TGraph*) fin->Get("roc_1b_r3");
  auto *g2b = (TGraph*) fin->Get("roc_2b_r3");
  auto *g2bsv = (TGraph*) fin->Get("roc_2b2sv_r3");
  auto *gAll = (TGraph*) fin->Get("roc_all_r3");
  

  auto hFrame = new TH1F ("hFrame","hFrame",1,0.0,1);
  hFrame->GetXaxis()->SetRangeUser(0.,1.);
  hFrame->GetYaxis()->SetRangeUser(0.00005,1.);
  hFrame->GetXaxis()->SetTitle("b-tag efficiency");
  hFrame->GetYaxis()->SetTitle("mistag rate");
  hFrame->GetXaxis()->CenterTitle(1);
  hFrame->GetYaxis()->CenterTitle(1);
  
  TCanvas *c=new TCanvas("c","c",600,600);
  c->SetGridx();
  c->SetGridy();
  c->SetLogy();


  g1b->SetLineColor(kblue);
  g2b->SetLineColor(kred);
  g2bsv->SetLineColor(kRed-9);
  gAll->SetLineColor(kgreen);

  g1b->SetLineWidth(6);
  g2b->SetLineWidth(6);
  g2bsv->SetLineWidth(6);
  gAll->SetLineWidth(6);
  
  
  hFrame->Draw();
  g1b->Draw("l");
  g2b->Draw("l");
  g2bsv->Draw("l");
  gAll->Draw("l");

  c->RedrawAxis();
  
  double yline = 0.001;  // target mistag rate
  
  
  TLine *line = new TLine(0., yline, 1., yline);
  //line->SetLineColor(kRed);
  line->SetLineStyle(2);  // dashed
  line->SetLineWidth(2);
  line->Draw("same");

  
  
  TLatex latex;
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.05,
		  yline*1.05,
		  "UParT-only mistag rate");


  TLatex latex2;
  latex2.SetTextSize(0.04);
  latex2.DrawLatex(0.025, 0.5,"100 < p_{T,jet} < 120 GeV, |#eta_{jet}| < 2");
  // hard-code sig. eff. and eff. mistag rate
  
  TMarker *wp = new TMarker(0.560412,0.000250118, 20);                                                                                                                                                            
  wp->SetMarkerSize(1.2);                                                                                                                                                                                         
  wp->SetMarkerColor(kRed-9);
  wp->Draw("same"); 
  /*
  TMarker *wp2 = new TMarker(0.560412,0.000250118, 24);                                                                                                                                                           
  wp2->SetMarkerSize(1.2);                                                                                                                                                                                       
  wp2->SetMarkerColor(kRed);
  wp2->Draw("same"); 
  */

  
  TLegend* leg = new TLegend(0.2, 0.55, 0.45, 0.8);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("UParT");
  leg->AddEntry(gAll, "all B (1B + #geq2B)", "l");
  leg->AddEntry(g1b,   "1B", "l");
  leg->AddEntry(g2b,   "#geq2B", "l");
  leg->AddEntry(g2bsv,   "#geq2B + 2SV", "l");
  leg->AddEntry(wp,"Analysis working point","p");
  leg->Draw();


  
  
  drawCMSLabel(c, "Internal Simulation", "2024 pp (5.36 TeV)");
}
