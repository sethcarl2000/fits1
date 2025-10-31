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

constexpr double mean_draw_range[] = { -8., 8. }; 


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

//NORMALIZED gaussian dist
double normalized_gauss(double x, double mean, double sigma)
{
  return 0.398942280401 * exp( - 0.5 * pow((x - mean)/sigma, 2) ) / sigma; 
} 


void fit1c(int entries=1000, bool save=false) 
{
  TRandom3 rand;

  const double mean  = rand.Uniform( mean_range[0], mean_range[1] );
  const double sigma = rand.Uniform( sigma_range[0], sigma_range[1] ); 

  //first, we must create the histogram and fill it with data according to a gaussian dist. 
  auto h = TH1D("h", "Test of Gauss dist;x;counts", n_bins, range[0], range[1]); 
  for (int i=0; i<n_points; i++) h.Fill( rand.Gaus()*sigma + mean ); 

  //now, we will do a log-liklihood fit of this data. 
  auto fitptr = h.Fit("gaus", "Q S N L");

  //let's caclulate our own log-liklihood...

  // uiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii12
  // --comment from assistant developer, Muon 

  const double LL_fit = fitptr->Chi2(); 

  using namespace std::placeholders; 
  auto norm_gaus = std::bind( normalized_gauss, _1, mean, sigma ); 

  const double LL = -log_liklihood(&h, norm_gaus); 

  auto h_LL = new TH1D("h_LL", "Distribution of Log-liklihood; (negative) Log Liklihood;", n_bins, 100., 175.); 

  printf(
    "Log liklihood calculated for %i points and %i bins:\n"
    " - % .2e (from fit)\n"
    " - % .2e (from my log_liklihood fcn)\n",
    n_points, n_bins, 
    LL_fit, 
    LL
  ); 

  //now, repeat this experiment many times to compare this value for reasonability. 
  for (int i=0; i<n_experiments; i++) {
    
    h.Reset(); 

    for (int i=0; i<n_points; i++) h.Fill( rand.Gaus()*sigma + mean ); 

    h_LL->Fill( -log_liklihood(&h, norm_gaus) ); 
  }
  
  auto c = new TCanvas("c", "Log-Liklihood comparison", 1400, 600);
  c->Divide(2,1); 
  
  
  c->cd(1); 
  h_LL->Draw("E, HIST"); 

  auto line = new TLine(LL,0., LL,h_LL->GetMaximum()); 
  line->SetLineColor(kRed); 
  line->Draw(); 

  //draw the CDF of our LL distribution 
  vector<double> x_LL, y_CDF; 
  const double total_stats = h_LL->Integral(); 
  auto ax = h_LL->GetXaxis(); 
  double cum=0.; 
  for (int b=1; b<=ax->GetNbins(); b++) {
    cum += h_LL->GetBinContent(b)/total_stats; 
    x_LL .push_back( ax->GetBinCenter(b) ); 
    y_CDF.push_back( 1. - cum ); 
  }

  c->cd(2); 
  auto g_CDF = new TGraph( x_LL.size(), x_LL.data(), y_CDF.data() ); 
  g_CDF->SetLineColor(kBlack);
  g_CDF->SetTitle("P-value from LL dist;Log Liklihood (LL);P(LL)"); 
  g_CDF->Draw(); 

  line = new TLine(LL,0., LL,1.); 
  line->SetLineColor(kRed); 
  line->Draw(); 

  return; 
}