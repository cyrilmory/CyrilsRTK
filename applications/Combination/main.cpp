#include <vector>
#include <string>
#include <sstream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "rtkSpectralFilter.h"
#include "rtkGgoFunctions.h"
#include "spectral_ggo.h"

using namespace itk;
using namespace  std;

int main(int argc, char* argv[])
{
    GGO(spectral, args_info);
    rtk::SpectralFilter::Image3D::Pointer ret = rtk::SpectralFilter::Image3D::New();
    rtk::SpectralFilter::Image3D::RegionType region;
    rtk::SpectralFilter::Image3D::SizeType size;
    std::vector<rtk::SpectralFilter::Image3D::Pointer> images;
    images.resize(args_info.inputImage_given);

    std::cout<<args_info.outputName_arg<<std::endl<<args_info.inputImage_given<<std::endl;
    for(int i = 0; i< args_info.inputImage_given; i++)
        images[i] = rtk::SpectralFilter::OpenImage3D(args_info.inputImage_arg[i]);

    //set size
    size[0] = images[0]->GetLargestPossibleRegion().GetSize()[0];
    size[1] = images[0]->GetLargestPossibleRegion().GetSize()[1];
    size[2] = args_info.inputImage_given;

    region.SetSize(size);
    ret->SetRegions(region);
    ret->Allocate();


    itk::ImageRegionIteratorWithIndex<rtk::SpectralFilter::Image3D> it(images[0], images[0]->GetLargestPossibleRegion());

    rtk::SpectralFilter::Image3D::IndexType index;
    for(it ; !it.IsAtEnd(); ++it)
    {
        index = it.GetIndex();
        for(int i = 0; i < images.size(); i++)
        {
            double val = (*images[i])[it.GetIndex()];
            index[2] = i;
            (*ret)[index] = val;
        }
    }
    rtk::SpectralFilter::SaveImage3D(ret, std::string(args_info.outputName_arg) +".mhd");
    return 0;
}
