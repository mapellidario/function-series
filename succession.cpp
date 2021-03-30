/*
Fourier Series
TODO: rendere logaritmiche le gif (non aggiungere un'immagine ogni n step,
ma aggiungere n immagini per decade)

g++ -o four fourier-series.cpp `root-config --cflags --glibs`
./four

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

// root
// tt le librerie includerli in "",
// altrimenti la va a cercare tra le librerie di sistema,
// ma io le ho salvate in un altra cartella.
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
// #include "TThread.h"
// #include "TStopwatch.h"
// #include "TStyle.h"
#include "TCanvas.h"
// #include "TFrame.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"

// double fourierk (double *x, double *par);

double fourierk (double *x, double *par) {
	double pi = TMath::Pi();
	double k = par[0] ;

	double y = sin(k * x[0]) + 1.5;
	return y;
}

int main (int numArg, char * listArg[])
{
	// constants
	double pi = TMath::Pi();

	// Setting
	double a0_half = 0. ;
	double xmin = -pi; //estremi del dominio
	double xmax = pi;
	int points = 100000; //numero di punti in cui calcolo la funzione nel dominio
	int iter = 11; // numero di termini della serie di fourier
	int start = 1; // minimo valore di k
					// quindi ho: \sum_{k=start}^iter

	//Setting Graphics
	const bool BatchPrint = 1 ;
	const bool gif = 1 ; //gif=1 stampa la gif animata.
						// NB è necessario che BatchPrint sia 1
	int GraphStep = 2 ; //ogni quante iterazioni stampare il grafico/
						//aggiungere un'immagine all'animazione gif
	double ymin = -2;
	double ymax = 2;

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
		fx[i]= a0_half ;//a_0/2 va segnato qui  //SETTARE
	}
	// vero e proprio ciclo: k per le iterazioni su termini di fourier,
	// i per le iterazioni sul dominio
	for (int k=start; k<iter; k++) {
		f[k] = new TF1 ("f1", fourierk, xmin, xmax, 5);
		f[k]->SetParameter(0, k);
		for (int i=0; i<points; i++) {
			fx[i] = f[k]->Eval(x[i]) ;
		} // parallelizzabile
		if (k%GraphStep==0) {
			for (int i=0; i<=points; i++) {
//				std::cout << i << " " << x[i] << " " << fx[i] <<  std::endl;
				graph->SetPoint(i,x[i],fx[i]);
				graph->SetTitle(Form(
					"f_{k=%d}(x) ",
					k));
				// graph->GetYaxis()->SetRangeUser(ymin,ymax); //SLOOOOOOOW
				graph->SetMarkerStyle(7);
			}
			graph->Draw("APL");
			if (BatchPrint && gif) {
				// gif+N: delay N*10 ms tra un immagine e l'altra
				std::string title_gif_partial =  aux_str + ".gif+20" ;
				c1->Print(title_gif_partial.c_str());
				std::cout << title_gif_partial << std::endl ;
			}
		}
		std::cout << "k: " << k << std::endl ;
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

/*Archive:
Square:
	int start = 1; //minimo valore di k
	for (int k=start; k<=iter; k++) {
		double normalize = 4./pi; //settare qui coefficiente moltiplicativo
		double ak = 0; //Settare qui ak e bk
		double bk = (double) 1 / (2*k-1) ;
		double coeff_cos = 0;
		double coeff_sin = (2*k-1);

triangle
	int start = 0; //minimo valore di k	//SETTARE
	for (int k=start; k<=iter; k++) {	//SETTARE
		double normalize = 8./pi/pi; //SETTARE
		double ak = 0; //Settare qui ak e bk //SETTARE
		double bk = (double) pow(-1,k) / (2*k+1)/(2*k+1) ; //SETTARE
		double coeff_cos = 0; //SETTARE
		double coeff_sin = (2*k+1); //SETTARE


sawtooth
	int start = 1; //minimo valore di k
	for (int k=start; k<=iter; k++) {
		double normalize = 2./pi;
		double ak = 0;
		double bk = (double) pow(-1,k) / k ;
		double coeff_cos = 0;
		double coeff_sin = k;

esercizio felli:
	for (int i=0; i<=points; i++) {
		x[i] = (double)(xmax-xmin)*i/points + xmin ;
		fx[i]=pi/4;//a_0/2 va segnato qui  //SETTARE
	}
	// vero e proprio ciclo: k per le iterazioni su termini di fourier,
	// i per le iterazioni sul dominio
	int start = 1; //minimo valore di k  //SETTARE
	for (int k=start; k<=iter; k++) { //SETTARE
		double normalize = 1; //settare qui coefficiente moltiplicativo
		double ak = (double) (1+pow(-1,k+1))/(pi*k*k); //Settare qui ak e bk
		double bk = (double) (1+2*pow(-1,k+1) )/ k ;
		double coeff_cos = k; //SETTARE
		double coeff_sin = k; //SETTARE

esercizio lucie: f(x) = 0 if x \in [-pi, 0], = -cos(x) if x \in [0, pi]

int start = 1; //minimo valore di k
for (int k=start; k<=iter; k++) {
	double normalize = 1./pi; //settare qui coefficiente moltiplicativo
	double ak = 0; //Settare qui ak e bk
	double bk = (double) - ( (pow(-1,k)+1) * k / (k*k-1) ) ;
	if (k==1) ak = - 1.570795 ;
	if (k==1) bk = 0. ;
	double coeff_cos = 0;
	double coeff_sin = k ;
	if (k==1) coeff_cos = 1. ;
	if (k==1) coeff_sin = 0. ;

esame felli 16/09/2015: abs(x)

*/
