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
    std::vector<rtk::SpectralFilter::Image3D::Pointer> vm;
    std::vector<vnl_vector<double> > vep;

    for(int i = 0; i < args_info.inputImage_given; i++)
        vm.push_back(rtk::SpectralFilter::OpenImage3D(args_info.inputImage_arg[i]));

    for(int i = 0; i < args_info.inputThickness_given; i++)
        vep.push_back(rtk::SpectralFilter::OpenThickness(args_info.inputThickness_arg[i]));


    double deg = 2;
    if(args_info.degree_given != 0)
        deg = args_info.degree_arg;

    std::vector<rtk::SpectralFilter::Image3D::Pointer> voutput;

    for(int i = 0; i < vep.size(); i++)
        voutput.push_back(s.Decomposition(vep[i], vm, deg));

    for(int i = 0; i < voutput.size() ; i++)
    {
        std::ostringstream number;
        number<<i;

        rtk::SpectralFilter::SaveImage3D(voutput[i], std::string(args_info.outputDir_arg) + "/polynomial_output_"+number.str()+".mhd");
        std::cout<<std::string(args_info.outputDir_arg) + "/polynomial_output_"+number.str()+".mhd"<<std::endl;
        //s.SeeImageValue(voutput[i]);
    }

    //cout<<"Sizeof = "<<sizeof(double)<<endl;

    return 0;
}
