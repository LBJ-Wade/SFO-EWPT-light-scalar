#include <iostream>
#include "simplebounce.h"
using namespace std;
using namespace simplebounce;

class MyModel : public GenericModel{
    public:
        MyModel(){
            setNphi(1);
        }
        double vpot (const double* phi) const{
            return phi[0]*phi[0]/2. - phi[0]*phi[0]*phi[0]*phi[0]/4.;
        }
        void calcDvdphi(const double *phi, double *dvdphi) const{
            dvdphi[0] = phi[0] - phi[0]*phi[0]*phi[0];
        }
};

int main() {
    BounceCalculator bounce;
    bounce.verboseOn();
    bounce.setRmax(1.);
    bounce.setDimension(4);
    bounce.setN(100);
    MyModel model;
    bounce.setModel(&model);

    double phiTV[1] = {10.};
    double phiFV[1] = {1.};
    bounce.setVacuum(phiTV, phiTV);

    bounce.solve();
    bounce.printBounce();

    cout << "S_E = " << bounce.action() << endl;

    return 0;
}
