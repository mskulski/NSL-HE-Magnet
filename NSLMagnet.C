#include "TCanvas.h"
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TF1.h>
#include "TH1F.h"
#include "TStyle.h"
#include "TColor.h"
#include <sstream> 
#include <fstream>
#include <stdlib.h>
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
using namespace std;

static int n = 10000;
Int_t steps = 1000;

Double_t am = 25*TMath::Pi()/180;
Double_t c = 299792458;
Double_t mu = 931.494*1e6/(c*c);
TRandom3* r3 = new TRandom3();

Double_t phi0;
Double_t theta0;
Double_t x0;
Double_t y0;

vector<vector<Double_t>> z(6);
vector<vector<Double_t>> x(6);
vector<vector<Double_t>> zo(6);
vector<vector<Double_t>> zf(6);
vector<vector<Double_t>> V(6);
vector<Double_t> m;
vector<Double_t> abun;

void Isotopes()
{
	NSLMagnet(3.5e6,3,3,"C",14,0,"");
	NSLMagnet(8.2e6,8,8,"Cl",36,0,"");
	NSLMagnet(7e6,9,9,"I",129,0,"");
	NSLMagnet(9.5e6,14,14,"Zr",93,16,"O");
	NSLMagnet(8.5e6,9,16,"Fe",60,16,"O");
}
	
void NSLMagnet(Double_t TV,Double_t q1,Double_t q,TString el,Double_t m0,Double_t m_add,TString el_add)
{	
	ifstream itable;
	itable.open("IsotopeTable.txt");
	Double_t icount;
	Double_t per;
	TString name1;
	Double_t massi;
	Double_t A1;
	Double_t per1;
	Double_t m0in;
	vector<TString> names;
	Int_t ntot = 0;
	Double_t amin = 1.01;
	Double_t amax = 0;
	Int_t m0max;
	
	TH1F *mass1[6];
	TH1F *mass12[6];
	TH1F *mass123[6];
	TGraph* graph[6];
	TMultiGraph *mg = new TMultiGraph();
	
	auto legend = new TLegend(0.1,0.75,0.4,0.9);
	auto legend1 = new TLegend(0.1,0.75,0.4,0.9);
	
	stringstream imgtitle;
	imgtitle << m0 << el;
	TString m_im = imgtitle.str();
	
	while(1)
	{
		itable >> name1 >> A1 >> massi >> per1;
		if(name1 == el)
		{	
			if(per1!=0 || A1==m0)
			{
				m.push_back(massi*mu);
				abun.push_back(per1);
				names.push_back(name1);
				if(A1==m0)
				{
					m0=massi*mu;
					m0in=icount;
				}
				icount++;
			}	
		}
		if(name1==el_add && A1==m_add)
		{
			m_add=massi*mu;
		}
		if(!itable.good())
		{
			break;
		}
	}
	itable.close();
	
	for(int s=0;s<icount;s++)
	{	
		if(s==m0in)
		{
			continue;
		}
		if(abun[s]<amin)
		{
			amin=abun[s];
		}
		if(abun[s]>amax)
		{	
			amax=abun[s];
			m0max=s;
		}
	}
	
	Double_t e0 = TV*((m0/(m0+m_add))+(q1/2)+(q/2));
	Double_t mp=m0;
	
	for(int j=0;j<n;j++)
	{
		y0 = r3->Gaus(0,0.002);
		x0 = r3->Gaus(0,0.002);
		if(ab(y0)>0.005 || ab(x0)>0.00075){continue;}
		phi0 = r3->Gaus(pi()/2,0.001);
		theta0 = r3->Gaus(pi()/2,0.001);
		Int_t vcount = 0;

		for(int s=0;s<icount;s++)
		{
			mp = m[s];
			Double_t e = r3->Gaus(TV*((mp/(mp+m_add))+(q1/2)+(q/2)),((mp/(mp+m_add))+(q1/2)+(q/2))*5e3);
			Double_t v0 = sqrt(2*e/mp);
			Double_t B = (1/1.016)*sqrt(2*m0*e0/(q*q));
			Double_t f = q*B/mp;
			zcalc(mp,v0,f,q,s,j,m0);
		}
		ntot++;
	}

	Double_t min = *min_element(zf[0].begin(),zf[0].end());
	Double_t max = *max_element(zf[icount-1].begin(),zf[icount-1].end());
	Double_t nbins = (max-min)/0.0001 + 1;

	TCanvas* mass = new TCanvas("Scan","Scan",3600,2700);
	
	for(int s=0;s<icount;s++)
	{
		mass1[s]=new TH1F(names[s],"",nbins,min,max);
		mass12[s]=new TH1F(names[s],"",nbins,min,max);
		
		for(int k=0;k<zf[s].size();k++)
		{	
			mass1[s]->Fill(zo[s][k]);
			mass12[s]->Fill(zf[s][k]);
		}
		
		Double_t amult = abun[s]/amin;
		mass1[s]->Scale(amult);
		mass12[s]->Scale(amult);
		
		stringstream namel;
		Int_t mround = TMath::Nint(m[s]/mu);
		namel << mround << names[s] << ":    " << mass1[s]->GetMean() << ",    "  << mass12[s]->GetMean();
		TString namel1 = namel.str();
		legend->AddEntry(mass1[s],namel1,"f");
		
		mass1[s]->SetFillColorAlpha(18+5*s,0.8);
		mass12[s]->SetFillColorAlpha(18+5*s,0.8);
		
		if(s==0)
		{
			mass12[0]->SetStats(kFALSE);
			mass12[0]->GetXaxis()->SetTitle("Magnet Exit/Max Position (m)");
			mass12[0]->GetXaxis()->CenterTitle();
			mass12[0]->GetXaxis()->SetTitleOffset(1.4);
			mass12[0]->GetXaxis()->SetLabelSize(.025);
			mass12[0]->SetTitle(m_im+" Magnet Exit vs Max FC Placement");
			mass12[0]->Draw("");
		}
		else{mass12[s]->Draw("SAME");}
		mass1[s]->Draw("SAME");
	}
	legend->Draw();

	Double_t maxy = mass1[m0max]->GetMaximum();
	
	for(int s=0;s<icount;s++)
	{
		if(maxy!=0)
		{
			mass1[s]->Scale(1/maxy);
			mass12[s]->Scale(1/maxy);
		}
		mass1[s]->SetAxisRange(0,1.25,"Y");
		mass12[s]->SetAxisRange(0,1.25,"Y");
	}
	
	Double_t min2 = 1e6;
	Double_t max2 = -1e6;
	
	for(int s=0;s<icount;s++)
	{
		Double_t mint = *min_element(V[s].begin(),V[s].end());
		Double_t maxt = *max_element(V[s].begin(),V[s].end());
		if(mint < min2){min2 = mint;}
		if(maxt > max2){max2 = maxt;}
	}
	Double_t nbins2 = abs(max2-min2)+1;
	
	TCanvas* mass2 = new TCanvas("Scan","Scan",3600,2700);
	
	for(int s=0;s<icount;s++)
	{
		mass123[s]=new TH1F(names[s],"",nbins2,min2,max2);
		mass123[s]->SetFillColorAlpha(18+5*s,0.8);
		
		for(int k=0;k<V[s].size();k++)
		{	
			if(s!=m0in)
			{
				mass123[s]->Fill(V[s][k]);
			}
		}
		
		stringstream namel;
		Int_t mround = TMath::Nint(m[s]/mu);
		namel << mround << names[s] << ":    " << mass123[s]->GetMean();
		TString namel1 = namel.str();
		legend1->AddEntry(mass123[s],namel1,"f");
		
		if(s==0)
		{
			mass123[0]->SetStats(kFALSE);
			mass123[0]->GetXaxis()->SetTitle("Terminal Voltage Difference Putting Beam at +/- 5.08 cm at Max Range (kV)");
			mass123[0]->GetXaxis()->CenterTitle();
			mass123[0]->GetXaxis()->SetTitleOffset(1.4);
			mass123[0]->GetXaxis()->SetLabelSize(.025);
			mass123[0]->SetTitle(m_im+" Terminal Voltage Differences for Stable Isotopes");
			mass123[0]->Draw("");
		}
		else{mass123[s]->Draw("SAME");}
	}
	legend1->Draw();	
	Double_t maxy1 = mass123[m0max]->GetMaximum();
	
	for(int s=0;s<icount;s++)
	{
		mass123[s]->SetAxisRange(0,1.25*maxy1,"Y");
	}
	
	TCanvas* mass3 = new TCanvas("Scan","Scan",3600,2700);
	
	for(int s=0;s<icount;s++)
	{
		graph[s] = new TGraph(ntot*steps*2,&x[s][0],&z[s][0]);
		graph[s]->SetMarkerColor(18+5*s);
		mg->Add(graph[s]);		
	}
	
	mg->SetTitle(m_im+" and Stable Isotope Trajectories");
	mg->Draw("AP");
	mg->GetXaxis()->SetLimits(-1.35,0.85);
	mg->SetMaximum(0.4);
	mg->SetMinimum(-1.25);
	
	TF1* mag1 = new TF1("mag1","-1.016 - tan(25*TMath::Pi()/180)*(1+x)",-1.266,-0.766);
	mag1->SetLineColor(1);
	mag1->Draw("same");
	TF1* mag2 = new TF1("mag2","-x/tan(25*TMath::Pi()/180)",-tan(am)*0.25,tan(am)*0.25);
	mag2->SetLineColor(1);
	mag2->Draw("same");
	TF1* mags1 = new TF1("mags1","0.0508",0.725,0.775);
	mags1->Draw("same");
	TF1* mags2 = new TF1("mags2","-0.0508",0.725,0.775);
	mags2->Draw("same");
	
	TImage *img = TImage::Create();
	img->FromPad(mass);
	stringstream imgname;
	imgname << m_im << ".png";
	TString imgf = imgname.str();
	img->WriteImage(imgf);
	
	TImage *img2 = TImage::Create();
	img2->FromPad(mass2);
	stringstream imgname2;
	imgname2 << m_im << "_V.png";
	TString imgf2 = imgname2.str();
	img2->WriteImage(imgf2);

	TImage *img3 = TImage::Create();
	img3->FromPad(mass3);
	stringstream imgname3;
	imgname3 << m_im << "_Trajectories.png";
	TString imgf3 = imgname3.str();
	img3->WriteImage(imgf3);
	
	for(int s=0;s<icount;s++)
	{
		zo[s].clear();
		zf[s].clear();
		V[s].clear();
		x[s].clear();
		z[s].clear();
	}
	m.clear();
	abun.clear();
}

void zcalc(Double_t mc,Double_t vc,Double_t fc,Double_t qi,Int_t v,Int_t j1,Double_t m0i)
{
	Double_t mc = m[v];
	Double_t K = 0.5*mc*vc*vc;
	Double_t Vc = 0;
	Double_t Vold = 1e5;
	Double_t Vmin = -1e6;
	Double_t msign = (mc-m0i)/ab(mc-m0i);
	Double_t Vmax = 1e6;
	Int_t cnt = 0;
	Int_t breakc = 0;
	
	while(breakc==0)
	{
	vc = sqrt(2*K/mc)*sqrt(1 + qi*Vc/K);
	Double_t vp = ab(vc*sn(theta0));
	Double_t vy = vc*cs(theta0);
	Double_t vx0 = vc*sn(theta0)*cs(phi0);
	Double_t t0 = -asn(vx0/vp)/fc;
	Double_t L = 2*1.016;

	Double_t p1 = vp/fc;
	Double_t p2 = tn(am);
 	Double_t p31 = 3.048 - 1.016*(1+tn(am)+tn(phi0)) + x0*tn(phi0);
  	Double_t p32 = p31/(tn(phi0)+tn(am));
	Double_t p3 = p32 + p1*cs(-fc*t0);
	Double_t p4 = -p32*tn(am) - 1.016*(1+tn(am)) - p1*sn(-fc*t0);
	Double_t oo = 2*(atn((sqrt(p1*p1*p2*p2 + p1*p1 - p2*p2*p4*p4 - 2*p2*p3*p4 - p3*p3) - p1*p2)/(p1 + p2*p4 + p3)));
	
	Double_t vzf = vp*cs(oo);
	Double_t vxf = vp*sn(oo);
	Double_t phif = atn(vzf/vxf);

	if(Vc==0)
	{
	Int_t stepstot = steps+1;
	for(int s=0;s<stepstot;s++)
	{
		z[v].push_back(p4 + p1*sn(s*oo/steps));
		x[v].push_back(p3 - p1*cs(s*oo/steps));
	}
	}
	
	Double_t x2 = p3 - p1*cs(oo);
	Double_t z2 = p4 + p1*sn(oo) + (0.0508*tn(am)-x2)*tn(phif);
	
	if(Vc==0){zo[v].push_back(z2);}
	
	Double_t Lf = 29.5*0.0254;
	Double_t t = Lf/vxf;
	Double_t zmfc = z2 + vzf*t;
	
	if(Vc==0)
	{
	Int_t stepstot = steps+1;
	for(int s=0;s<stepstot;s++)
	{
		z[v].push_back(z2 + vzf*t*s/steps);
		x[v].push_back(vxf*t*s/steps);
	}
	}
	
	if(Vc==0){zf[v].push_back(zmfc);}
	
	if(Vc == Vold){cnt++;}
	else{Vold = Vc;}
	
	if(ab(zmfc-msign*0.0508) <= 1e-12 || cnt > 10)
	{
		V[v].push_back(Vc/1000);
		breakc = 1;
	}
	else if(zmfc-msign*0.0508 > 0)
	{
		Vmax = Vc;
		Vc = (Vmax + Vmin)/2;
	}
	else if(zmfc-msign*0.0508 < 0)
	{
		Vmin = Vc;
		Vc = (Vmax + Vmin)/2;
	}
	}
}
Double_t tn(double angle)
{
	return TMath::Tan(angle);
}
Double_t atn(double angle)
{
	return TMath::ATan(angle);
}
Double_t asn(double angle)
{
	return TMath::ASin(angle);
}
Double_t sn(double angle)
{
	return TMath::Sin(angle);
}
Double_t cs(double angle)
{
	return TMath::Cos(angle);
}
Double_t acs(double angle)
{
	return TMath::ACos(angle);
}
Double_t ab(double xv)
{
	return TMath::Abs(xv);
}
Double_t pi()
{
	return TMath::Pi();
}
