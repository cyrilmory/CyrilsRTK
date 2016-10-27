#include <vector>
#include <string>
#include <sstream>
#include "rtkSpectralFilter.h"
#include "rtkGgoFunctions.h"
#include "spectral_ggo.h"


int main(int argc, char* argv[])
{
    GGO(spectral, args_info);

    rtk::SpectralFilter s;
    std::vector<rtk::SpectralFilter::Image2D::Pointer> vm;
    std::vector<rtk::SpectralFilter::Image3D::Pointer> vp;
    std::vector<rtk::SpectralFilter::Image2D::Pointer> voutput;
    for(int i = 0; i < args_info.inputImage_given; i++)
        vm.push_back(rtk::SpectralFilter::OpenImage2D(args_info.inputImage_arg[i]));

    for(int i = 0; i < args_info.inputPolynomial_given; i++)
        vp.push_back(rtk::SpectralFilter::OpenImage3D(args_info.inputPolynomial_arg[i]));

    double deg = 2;
    if(args_info.degree_given != 0)
        deg = args_info.degree_arg;

    for(int i = 0; i < vp.size(); i++)
        voutput.push_back(s.Forward(vp[i], vm, deg));

    for(int i = 0; i < voutput.size(); i++)
    {
        std::ostringstream number;
        number<<i;
        rtk::SpectralFilter::SaveImage2D(voutput[i], std::string(args_info.outputDir_arg) + "/image_output_"+number.str()+".mhd");
        std::cout<<std::string(args_info.outputDir_arg) + "/image_output_"+number.str()+".mhd"<<std::endl;
    }
    return 0;
}
