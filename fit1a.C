#include <TH1D.h> 
#include <TCanvas.h> 
#include <TPad.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
#include <TRandom3.h> 
#include <TScatter.h> 
#include <vector> 
#include <TStyle.h> 
#include <cstdio> 
#include <cmath> 
#include <iostream> 

using namespace std; 

const int n_points = 1e3; 
const int n_experiments = 1e4; 
const int n_bins = 100; 
constexpr double range[] = {0., 100.}; 

// the maximum and minimum mean and sigma of any 'experiment'. 
// the mean and sigma of each experiment will be generated uniformly over these ranges 
const double mean_range[] = { 50., 50. };
const double sigma_range[] = { 10., 10. }; 

void fit1a(int entries=1000, bool save=false) 
{
  TRandom3 rand;

  auto h = TH1D("h", "Test of Gauss dist;x;counts", n_bins, range[0], range[1]); 
  
  auto h_mean = new TH1D("h_mean", "Mean (from fit) - mean (actual);#delta mean (fit - actual);", n_bins, -5., 5.); 
  auto h_chi2 = new TH1D("h_chi2", "Reduced #chi^{2} (from fit);#chi^{2} / N. DoF;", n_bins, 0., 2.5);
  auto h_merr = new TH1D("h_prob", "Error of mean, as reported by fit;#sigma(#bar{x});", n_bins, 0., 1.);  

  vector<double> reduced_chi2, p_chi2;  
  
  printf("conducting %i experiments, with %i pts each...", n_experiments, n_points); cout << flush; 
  for (int i=0; i<n_experiments; i++) {

    double mean = rand.Uniform(mean_range[0], mean_range[1]);
    double sigma = rand.Uniform(sigma_range[0], sigma_range[1]);  

    //fill the histogram
    for (int i=0; i<n_points; i++) h.Fill( mean + rand.Gaus()*sigma ); 

    //now, perform the fit
    auto fitptr = h.Fit("gaus", "Q S N");
    
    if (!fitptr->IsValid()) continue; 

    h_mean->Fill( fitptr->Parameter(1) - mean ); 

    h_chi2->Fill( fitptr->Chi2()/fitptr->Ndf() ); 
    h_merr->Fill( fitptr->ParError(1) ); 

    reduced_chi2.push_back( fitptr->Chi2()/fitptr->Ndf() ); 
    p_chi2      .push_back( fitptr->Prob() ); 

    h.Reset(); 
  }
  cout << "done." << endl; 

  auto c = new TCanvas("c", "result of fits", 1400, 800); 
  c->Divide(2,2); 

  auto h_prob_x2 = new TScatter(p_chi2.size(), reduced_chi2.data(), p_chi2.data(), nullptr, nullptr); 
  h_prob_x2->SetTitle("Reduced #chi^{2} vs P(#chi^{2});#chi^{2} / N.DoF;P(#chi^{2})");

  c->cd(1); h_chi2->SetStats(0); h_chi2->Draw("E"); 
  c->cd(2); h_mean->SetStats(1); h_mean->Draw("E, HIST");
  
  c->cd(3); 
  //h_prob_x2->GettYaxis()->SetRangeUser(0., h_prob->GetMaximum()*1.2); 
  //find min and max elements of a std::vector<double> 
  auto vec_min = [](const vector<double>& v) { 
    double min = +std::numeric_limits<double>::infinity(); 
    for (const auto& x : v) if (x < min) min = x; 
    return min; 
  };
  auto vec_max = [](const vector<double>& v) { 
    double max = -std::numeric_limits<double>::infinity(); 
    for (const auto& x : v) if (x > max) max = x; 
    return max; 
  };

  //h_prob_x2->SetMarkerStyle(kOpenCircle);
  h_prob_x2->SetMarkerSize(0.75); 
  h_prob_x2->Draw("SAME"); 

  c->cd(4); h_merr->SetStats(1); h_merr->Draw("E, HIST"); 

  return; 
}