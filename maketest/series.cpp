/*
Function Series
TODO: rendere logaritmiche le gif (non aggiungere un'immagine ogni n step,
ma aggiungere n immagini per decade)

ACHTUNG: TMath::Gamma()

g++ -o series series.cpp `root-config --cflags --glibs`
./series
*/

//std
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
//#include <ctime>
//#include <vector>

//root
//tt le librerie includerli in "", altrimenti la va a cercare tra le librerie di sistema, ma io le ho salvate in un altra cartella.
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TThread.h"
#include "TStopwatch.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TSlider.h"

#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TF2.h"

#include "TList.h"
#include "TArrayD.h"
#include "TString.h"

#include "TMath.h"
#include "TRandom.h"
#include "TVirtualFFT.h"

double functionk (double *x, double *par);
double factorial (double x) ;

//SETTARE: termine n-esimo della serie
double functionk (double *x, double *par) {
//	return (1. - fabs(x[0])) / (par[0]*par[0]) ; //x[i]
// return (1. - sqrt(x[0]*x[0]+1./(2.*pow(par[0],5)))) / (par[0]*par[0]) ; //x[i]
//    return (pow(-1,par[0]-1) / (par[0]) * pow(x[0] - 3, par[0] )) ;
    // test functions
    // x / (1 - x), x \in (-1,1)
//    return pow(x[0], par[0]) ;
    // e^x \forall x
//    return pow(x[0], par[0]) / factorial (par[0]) ; // OK
//    return pow(x[0], par[0]) / TMath::Factorial (par[0]) ; //OK
    // sin(x) \forall x
//    return pow (-1, par[0]) * pow (x[0], 2*par[0]+1) / TMath::Factorial(2*par[0]+1) ; //OK
}

// double factorial (double x) {
//  if ( x == 1) return 1 ;
//  else return x * factorial (x-1) ;
// }

int main (int numArg, char * listArg[])
{
	double pi = TMath::Pi();
	//Setting Graphics
	const bool BatchPrint = 0;
	const bool gif = 0; //gif=1 stampa la gif animata. NB è necessario che BatchPrint sia 1
	int points = 100; //numero di punti in cui calcolo la funzione nel dominio
	int iter = 150; //numerodi termini della serie di fourier
	int start = 0; //minimo valore di k  //SETTARE
	int GraphStep = 1; //ogni quante iterazioni stampare il grafico/aggiungere un'immagine all'animazione gif
	int GraphStepVec[10] = {0,1,2,3,4,5,6,7,8,9};
	int LogCoeff = 10;
	double xmin = -10; //estremi del dominio
	double xmax = 5 ;
	double ymin = -4.; //usati per il rangeuser sul graph
	double ymax = 4.;

//------------------------------------------------------------------------------

	//general declarations
	std::string aux_str = listArg[0] ;
	TApplication* theApp = new TApplication("App", &numArg, listArg);
	if (BatchPrint) {
		gROOT->SetBatch(); //Activate batch (edit images without opening them)
		if (gif) {
			std::string title_old_gif =  aux_str + ".gif" ;
			gSystem->Unlink(title_old_gif.c_str()); //dovrebbe essere title_gif
		}
	}
	TCanvas *c1 = new TCanvas ("c1","First Canvas",1900,1000);

	//generating the Graphs
	TGraph *graph = new TGraph();
	double x[points+1];
	double fx[points+1];
	TF1 *f[iter];
	//riempio il vettore delle x (i punti del dominio in cui valutare la funzione)
	//inizializzo fx per non avere brutte sorprese
	for (int i=0; i<=points; i++) {
		x[i] = (double)(xmax-xmin)*i/points + xmin ;
		fx[i]=0.; //SETTARE: offset
	}
	//vero e proprio ciclo: k per le iterazioni su termini di fourier, i per le iterazioni sul dominio
	for (int k=start; k<=iter; k++) { //SETTARE
		//double normalize = 1.; //settare qui coefficiente moltiplicativo //SETTARE
		std::cout << "k: " << k << std::endl ;
		f[k] = new TF1 ("f1", functionk, xmin, xmax, 1);
		f[k]->SetParameter(0,k);
//		std::cout << f[k]->GetParameter(0) << std::endl ;
		for (int i=0; i<=points; i++) {
			fx[i] += f[k]->Eval(x[i]) ;
            std::cout << "i = " << i << " x[i] = " << x[i] << " f(x) = " << fx[i] << std::endl ;;
		}
		if (k%GraphStep==0) {
//		if (k==GraphStepVec[1] || k==GraphStepVec[3] || k==GraphStepVec[5] ||
//										k==GraphStepVec[7] || k==GraphStepVec[9] ) {

			for (int i=0; i<=points; i++) {
//				std::cout << i << " " << x[i] << " " << fx[i] <<  std::endl;
				graph->SetPoint(i,x[i],fx[i]);
			}
			graph->SetTitle(Form("f(x)=#Sigma_{k=%d}^{%d}f_{k} ",start,k));
//			graph->GetYaxis()->SetRangeUser(ymin,ymax);
			graph->SetMarkerStyle(7);
			graph->Draw("APL");
			if (BatchPrint && gif) {
				std::string title_gif_partial =  aux_str + ".gif+20" ;
					//gif+N: delay N*10 ms tra un immagine e l'altra
				c1->Print(title_gif_partial.c_str());
				std::cout << title_gif_partial << std::endl ;
			}
			if (k==GraphStepVec[9]) {
				for (int j=0; j<10; j++){
					GraphStepVec[j] *= LogCoeff ;
					std::cout << GraphStepVec[j] << std::endl ;
				}
			}
		}
	}
//	TCanvas *c2 = new TCanvas ("c2","Second Canvas",1900,1000);
//	f[3]->SetNpx(1000);
//	f[3]->Draw();

//------------------------------------------------------------------------------
	if (BatchPrint) {
		//Stampa Canvas
		std::string aux = listArg[0] ;
		c1->Modified();
		c1->Update();
		c1->Draw();
		//decommentare per generare il file .png
		std::string title_png;
		std::string png = ".png";
		title_png = aux + png;
		c1->Print(title_png.c_str(),"png");
		//generare il file .pdf
		std::string title_pdf;
		std::string pdf = ".pdf";
		title_pdf = aux + pdf;
		c1->Print(title_pdf.c_str(),"pdf");
		//generare il file .tex
		std::string title_tex;
		std::string tex = ".tex";
		title_tex = aux + tex ;
		c1->Print(title_tex.c_str());
		//rendere infinito il loop sul file .gif
		if (gif) {
			std::string title_gif =  aux + ".gif++" ;
			c1->Print(title_gif.c_str());
		}
	}
	if (!BatchPrint) theApp->Run(); //Deactivate the application

	return 0;
}
