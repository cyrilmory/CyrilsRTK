#include <vector>
#include <string>
#include <sstream>
#include "rtkSpectralFilter.h"
#include "rtkGgoFunctions.h"
#include "spectral_ggo.h"

using namespace  std;

int main(int argc, char* argv[])
{
    GGO(spectral, args_info);

    rtk::SpectralFilter s;
    std::vector<rtk::SpectralFilter::Image3D::Pointer> vm;
    vnl_vector<double> ep;
    rtk::SpectralFilter::Image3D::Pointer pol;
    double deg = 2;

    for(int i = 0; i < args_info.inputImage_given; i++)
    {
        vm.push_back(rtk::SpectralFilter::OpenImage3D(args_info.inputImage_arg[i]));
    }

    ep = rtk::SpectralFilter::OpenThickness(args_info.inputThickness_arg);
    pol = rtk::SpectralFilter::OpenImage3D(args_info.inputPolynomial_arg);

    if(args_info.degree_given != 0)
        deg = args_info.degree_arg;


    s.ShowEstimatedError(vm, pol, ep, deg);

    return 0;
}
