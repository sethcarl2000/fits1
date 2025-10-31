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

const int n_points = 1e2; 
const int n_experiments = 1e4; 
const int n_bins = 100; 
constexpr double range[] = {0., 100.}; 

// the maximum and minimum mean and sigma of any 'experiment'. 
// the mean and sigma of each experiment will be generated uniformly over these ranges 
const double mean_range[] = { 50., 50. };
const double sigma_range[] = { 10., 10. }; 

constexpr double mean_draw_range[] = { -8., 8. }; 

void fit1b(int entries=1000, bool save=false) 
{
  TRandom3 rand;

  auto h = TH1D("h", "Test of Gauss dist;x;counts", n_bins, range[0], range[1]); 
  
  auto h_chi2 = new TH1D("h_mean", Form("Mean (from #chi^{2} fit): %i points;#Delta mean (#chi^{2}-fit - actual);",n_points), 
    n_bins, mean_draw_range[0], mean_draw_range[1]);
   
  auto h_LL   = new TH1D("h_LL",   Form("Mean (from LL fit): %i points;#Delta mean (LL-fit - actual);",n_points), 
    n_bins, mean_draw_range[0], mean_draw_range[1]); 
  
  printf("conducting %i experiments, with %i pts each...", n_experiments, n_points); cout << flush; 
  for (int i=0; i<n_experiments; i++) {

    double mean = rand.Uniform(mean_range[0], mean_range[1]);
    double sigma = rand.Uniform(sigma_range[0], sigma_range[1]);  

    //fill the histogram
    for (int i=0; i<n_points; i++) h.Fill( mean + rand.Gaus()*sigma ); 

    //now, perform the fit
    auto fp_chi2 = h.Fit("gaus", "Q S N");
    auto fp_LL   = h.Fit("gaus", "Q S N L");
    
    if (fp_chi2->IsValid()) h_chi2->Fill( fp_chi2->Parameter(1) - mean );   
    if (fp_LL  ->IsValid()) h_LL  ->Fill( fp_LL  ->Parameter(1) - mean ); 
 
    h.Reset(); 
  }
  cout << "done." << endl; 

  auto c = new TCanvas("c", "result of fits", 1400, 800); 
  c->Divide(2,1); 

  c->cd(1); h_chi2->SetStats(1); h_chi2->Draw("E"); 
  c->cd(2); h_LL  ->SetStats(1); h_LL  ->Draw("E");
  
  return; 
}