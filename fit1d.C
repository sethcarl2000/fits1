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
#include <functional> 
#include <iostream> 
#include <TAxis.h> 
#include <TLine.h> 
#include <TGraph.h> 

using namespace std; 

const int n_points = 1e2; 
const int n_experiments = 1e5; 
const int n_bins = 100; 
constexpr double range[] = {0., 100.}; 

// the maximum and minimum mean and sigma of any 'experiment'. 
// the mean and sigma of each experiment will be generated uniformly over these ranges 
const double mean_range[] = { 50., 50. };
const double sigma_range[] = { 10., 10. }; 

constexpr int    mean_draw_pts = 100.; 
constexpr double mean_draw_range = 8.; 

//assumes PDF is properly normailized! 
double log_liklihood(TH1D* hist, function<double(double)> PDF) 
{
  double LL = 0.; 

  //actual log of a factorial. let's see how slow this is... 
  auto logfact = [](int n) {
    double lf = 0.; if (n < 2) return lf;  
    for (int i=1; i<=n; i++) lf += log(n); 
    return lf; 
  };
  
  auto ax = hist->GetXaxis(); 
  
  const int nbins = ax->GetNbins(); 

  const double total_stats = hist->Integral(); 
  const double dx = ( ax->GetXmax() - ax->GetXmin() )/((double)nbins-1); 

  for (int b=1; b<=nbins; b++) {

    //center of this bin 
    double x = ax->GetBinCenter(b); 

    //average number of stats we should EXPECT in this bin (if PDF is correct)
    double expect = PDF(x) * dx * total_stats; 

    double bin_stats = hist->GetBinContent(b); 

    //now, add the log of the poisson prob. that 'N' stats would end up in this bin. 
    
    
    LL += (bin_stats * log(expect)) - expect - logfact(bin_stats); 
  }

  //return the total log liklihood
  return LL; 
} 

//assumes PDF is properly normailized! 
double chi2(TH1D* hist, function<double(double)> PDF) 
{
  double chi2 = 0.; 

  auto ax = hist->GetXaxis(); 
  
  const int nbins = ax->GetNbins(); 

  const double total_stats = hist->Integral(); 
  const double dx = ( ax->GetXmax() - ax->GetXmin() )/((double)nbins-1); 

  for (int b=1; b<=nbins; b++) {

    //center of this bin 
    double x = ax->GetBinCenter(b); 

    //average number of stats we should EXPECT in this bin (if PDF is correct)
    double expect = PDF(x) * dx * total_stats; 

    double bin_stats = hist->GetBinContent(b); 

    if (bin_stats < 0.1) continue;  
    //now, add the log of the poisson prob. that 'N' stats would end up in this bin. 
    
    chi2 += pow( (bin_stats - expect), 2 )/expect;  
  }

  //return the total log liklihood
  return chi2; 
} 

//NORMALIZED gaussian dist
double normalized_gauss(double x, double mean, double sigma)
{
  return 0.398942280401 * exp( - 0.5 * pow((x - mean)/sigma, 2) ) / sigma; 
} 


void fit1d(int entries=1000, bool save=false) 
{
  TRandom3 rand;

  const double mean  = rand.Uniform( mean_range[0], mean_range[1] );
  const double sigma = rand.Uniform( sigma_range[0], sigma_range[1] ); 

  //first, we must create the histogram and fill it with data according to a gaussian dist. 
  auto h25 = TH1D("h25", "Test of Gauss dist;x;counts", n_bins, range[0], range[1]); 
  auto h1k = TH1D("h1k", "Test of Gauss dist;x;counts", n_bins, range[0], range[1]); 
  
  TFitResultPtr fitptr_25, fitptr_1k;

  //repeat this until we get a valid fit (sometimes it might fail!)
  cout << "Trying to get valid fit..." << flush; 
  do {

    h25.Reset(); 
    h1k.Reset(); 

    //fill our histograms 
    for (int i=0; i<25; i++)  h25.Fill( rand.Gaus()*sigma + mean ); 
    for (int i=0; i<1e3; i++) h1k.Fill( rand.Gaus()*sigma + mean );

    fitptr_25 = h25.Fit("gaus", "S Q N L");
    fitptr_1k = h1k.Fit("gaus", "S Q N");  

  } while (!(fitptr_25->IsValid() && fitptr_1k->IsValid())); 
  cout << "done!" << endl; 

  //we're going to vary the mean slowly. 
  vector<double> pts_mean, pts_LL, pts_chi2; 

  double test_mean = mean - mean_draw_range; 

  const double dm = 2. * mean_draw_range / ((double)mean_draw_pts-1); 

  using namespace std::placeholders; 

  //now, repeat this experiment many times to compare this value for reasonability. 
  for (int i=0; i<mean_draw_pts; i++) {
    pts_mean.push_back( test_mean - mean ); 

    auto norm_gaus = std::bind(normalized_gauss, _1, test_mean, sigma); 

    pts_LL  .push_back( -log_liklihood(&h25, norm_gaus) ); 
    pts_chi2.push_back( chi2(&h1k, norm_gaus) ); 

    test_mean += dm; 
  }
  
  auto c = new TCanvas("c", "Dependence of LL and Chi2 on mean", 1400, 600);
  c->Divide(2,1); 
  
  auto vec_min = [](const vector<double>& v) {
    double m=+1.e100; for (auto x : v) if (x < m) m = x; return m;  
  }; 
  auto vec_max = [](const vector<double>& v) {
    double m=-1.e100; for (auto x : v) if (x > m) m = x; return m;  
  }; 

  TLine *line; 
  double fit_mean, fit_mean_err; 
  
  c->cd(1); 
  auto g_LL = new TGraph( pts_mean.size(), pts_mean.data(), pts_LL.data() ); 
  g_LL->SetTitle("Variance of -LL with mean (25 pts);Mean (diff. from actual);-Log Liklihood");
  g_LL->Draw(); 

  fit_mean      = fitptr_25->Parameter(1) - mean; 
  fit_mean_err  = fitptr_25->ParError(1); 
  line = new TLine(
    fit_mean, vec_min(pts_LL), 
    fit_mean, vec_max(pts_LL)
  ); 
  line->SetLineColor(kRed);
  line->Draw("SAME"); 
  line = new TLine(
    fit_mean + fit_mean_err, vec_min(pts_LL), 
    fit_mean + fit_mean_err, vec_max(pts_LL)
  ); 
  line->SetLineColor(kRed); line->SetLineStyle(kDotted); 
  line->Draw(); 
  line = new TLine(
    fit_mean - fit_mean_err, vec_min(pts_LL), 
    fit_mean - fit_mean_err, vec_max(pts_LL)
  );
  line->SetLineColor(kRed); line->SetLineStyle(kDotted); 
  line->Draw(); 


  c->cd(2); 
  auto g_chi2 = new TGraph( pts_mean.size(), pts_mean.data(), pts_chi2.data() ); 
  g_chi2->SetTitle("Variance of #chi^{2} with mean (10^{3} pts);Mean (diff. from actual);#chi^{2}");
  g_chi2->Draw(); 

  fit_mean      = fitptr_1k->Parameter(1) - mean; 
  fit_mean_err  = fitptr_1k->ParError(1); 
  line = new TLine(
    fit_mean, vec_min(pts_chi2), 
    fit_mean, vec_max(pts_chi2)
  ); 
  line->SetLineColor(kRed);
  line->Draw(); 
  line = new TLine(
    fit_mean + fit_mean_err, vec_min(pts_chi2), 
    fit_mean + fit_mean_err, vec_max(pts_chi2)
  ); 
  line->SetLineColor(kRed); line->SetLineStyle(kDotted); 
  line->Draw(); 
  line = new TLine(
    fit_mean - fit_mean_err, vec_min(pts_chi2), 
    fit_mean - fit_mean_err, vec_max(pts_chi2)
  );
  line->SetLineColor(kRed); line->SetLineStyle(kDotted); 
  line->Draw(); 

  return; 
}