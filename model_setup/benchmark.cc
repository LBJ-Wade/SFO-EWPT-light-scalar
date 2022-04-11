// table 1 (2 fields) of 1901.03714
#define _USE_MATH_DEFINES // for C++
#include<iostream>
#include<cmath>
#include"simplebounce.h"
#include<sys/time.h>
#include <numbers>
using namespace std;
using namespace simplebounce;

class MyModel : public GenericModel{
  public:
    double lambda;
    double muH;
    double muS;
    double A;
    double g;
    double gY;
    double yt;
    double Qsq;
	  MyModel(){
	  	  setNphi(2);
        lambda = 0.125299;
        A = .106204;
        muH = .202598;
        muS = .218138;
        g = 0.65;
        gY = 0.36;
        yt = 0.9945;
        Qsq = 1.50*1.50;
	}

	  // potential of scalar field(s)
	  double vpot(const double* phi) const{
      double vtree = (0.25*lambda*phi[0]*phi[0]*phi[0]*phi[0]-0.5*muH*muH*phi[0]*phi[0]+0.5*muS*muS*phi[1]*phi[1]-0.5*A*phi[1]*phi[0]*phi[0]);
      double mW = (g*g*phi[0]*phi[0]*0.25);
      double mZ = ((g*g+gY*gY)*phi[0]*phi[0]*0.25);
      double mt = (yt*yt*phi[0]*phi[0]*0.5);
      double vcw = (6*mW*mW*(log(mW/Qsq)-5./6)+3*mZ*mZ*(log(mZ/Qsq)-5./6)-12.*mt*mt*(log(mt/Qsq)-1.5))/(64*M_PI*M_PI);
	  	return vtree+vcw;
	}

	// derivative of potential of scalar field(s)
	  void calcDvdphi(const double* phi, double* dvdphi) const{
      double mW = (g*g*phi[0]*phi[0]*0.25);
      double mZ = ((g*g+gY*gY)*phi[0]*phi[0]*0.25);
      double mt = (yt*yt*phi[0]*phi[0]*0.5);
      double mWd = (g*g*phi[0]*0.5);
      double mZd = ((g*g+gY*gY)*phi[0]*0.5);
      double mtd = (yt*yt*phi[0]);
	  	dvdphi[0] = lambda*phi[0]*phi[0]*phi[0] - muH*muH*phi[0]-A*phi[1]*phi[0]+(6*mW*mWd*(2*log(mW/Qsq)-2./3)+3*mZ*mZd*(2*log(mZ/Qsq)-2./3)-12*mt*mtd*(2*log(mt/Qsq)-2.))/(64*M_PI*M_PI);
	  	dvdphi[1] = muS*muS*phi[1]-0.5*A*phi[0]*phi[0];
	  }

//    void dvdphi(const double* phi) const{
//      double dv[2];
//      double mW = (g*g*phi[0]*phi[0]*0.25);
//      double mZ = ((g*g+gY*gY)*phi[0]*phi[0]*0.25);
//      double mt = (yt*yt*phi[0]*phi[0]*0.5);
//      double mWd = (g*g*phi[0]*0.5);
//      double mZd = ((g*g+gY*gY)*phi[0]*0.5);
//      double mtd = (yt*yt*phi[0]);
//      dv[0] = lambda*phi[0]*phi[0]*phi[0] - muH*muH*phi[0]-A*phi[1]*phi[0]+(6*mW*mWd*(2*log(mW/Qsq)-2./3)+3*mZ*mZd*(2*log(mZ/Qsq)-2./3)-12*mt*mtd*(2*log(mt/Qsq)-2.))/(64*M_PI*M_PI);
//      dv[1] = muS*muS*phi[1]-0.5*A*phi[0]*phi[0];
//      cout << "dvdh = " << dv[0] << ", dvdS = " << dv[1] << endl;
//    }

};




int main() {

	BounceCalculator bounce;
  bounce.verboseOn();
	bounce.setRmax(1); // phi(rmax) = phi(False vacuum)
	bounce.setDimension(4); // number of space dimension
	bounce.setN(100); // number of grid
	MyModel model;
  //double location[2] = {86.,83.};
  //cout << model.vpot(location) << endl;
  //model.dvdphi(location);
	bounce.setModel(&model);

	double phiTV[2] = {15.,150.}; // a point at which V<0
	double phiFV[2] = {0.868666,0.838824}; // false vacuum
	bounce.setVacuum(phiTV, phiFV);

	// calcualte the bounce solution
	bounce.solve();

  bounce.printBounce();
  bounce.writeBounce("output/benchmark.csv");
	// Euclidean action
	cout << "S_E = " << bounce.action() << endl;

	return 0;
}
