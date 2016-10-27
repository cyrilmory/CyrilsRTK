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

    rtk::SpectralFilter::Image2D:: Pointer i1 = rtk::SpectralFilter::OpenImage2D(args_info.inputImage1_arg);
    rtk::SpectralFilter::Image2D:: Pointer i2 = rtk::SpectralFilter::OpenImage2D(args_info.inputImage2_arg);
    rtk::SpectralFilter::Image2D::Pointer ret = rtk::SpectralFilter::Image2D::New();
    ret->SetRegions(i1->GetLargestPossibleRegion());
    ret->Allocate();
    itk::ImageRegionIteratorWithIndex<rtk::SpectralFilter::Image2D> it(i1, i1->GetLargestPossibleRegion());

    for(;!it.IsAtEnd(); ++it)
    {
        (*ret)[it.GetIndex()] = (it.Get() - (*i2)[it.GetIndex()])/(*i2)[it.GetIndex()];
    }

    rtk::SpectralFilter::SaveImage2D(ret, args_info.outputName_arg);

    return 0;
}
