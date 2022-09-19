#include <iostream>
#include <cmath>
#include "WLP2D.h"

double Jz(double t){
    return  ((t-t0)/St)*exp(-(t - t0)*(t -t0) / (2.0*St*St)); 
}
