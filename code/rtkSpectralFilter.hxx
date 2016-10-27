#ifndef SPECTRAL_HXX
#define SPECTRAL_HXX

#include "rtkSpectralFilter.h"

namespace rtk
{

SpectralFilter::SpectralFilter()
{

}

//Public function

SpectralFilter::Image3D::Pointer SpectralFilter::Decomposition(vnl_vector<double>  ep,
                    const std::vector<Image3D::Pointer>& vm, double deg)
{
    Image3D::RegionType region1 = vm[0]->GetLargestPossibleRegion();
    Image3D::IndexType index;
    Image3D::SizeType size;

    size[0] = region1.GetSize()[0];
    size[1] = region1.GetSize()[1];
    /*size[0] = 361;
    size[1] = 167;*/

    //Return images initialization
    Image3D::Pointer ret = Image3D::New();
    Image3D::RegionType region;
    region.SetSize(size);
    size[2] = PolynomialSize(vm.size(), deg);
    region.SetSize(size);
    ret->SetRegions(region);
    ret->Allocate();

    itk::ImageRegionIteratorWithIndex<Image3D> it(vm[0], vm[0]->GetLargestPossibleRegion());

    for(; it.GetIndex()[2] == 0; ++it)
    {
        //Index Update
        index[0] = it.GetIndex()[0];
        index[1] = it.GetIndex()[1];

        //Construction of k pixel line (size P). with k number of energy level
        std::vector<vnl_vector<double> > vom;
        for(int k = 0; k < vm.size(); k++)
            vom.push_back(TakePixelVector(vm[k], index[0], index[1]));

        //Construction of am matrix for solving Ax=b with A = am
        //The am matrix is based on vom values
        vnl_sparse_matrix<double> am = ConstructMatrix(vom, deg);


        //Solve Ax=b for each energy level
        std::vector<double> x = SolveLinearSystem(am, ep);

        /*for(int i = 0; i < size[0]; i ++)
            for(int j = 0; j < size[1]; j++)
            {
                index[0] = i;
                index[1] = j;*/
                for(int l = 0; l < x.size(); l++)
                {
                    index[2] = l;
                    (*ret)[index] = x[l];
                }
//            }
    }

    return ret;

}

SpectralFilter::Image2D::Pointer SpectralFilter::Forward(Image3D::Pointer ip, std::vector<Image2D::Pointer> vm,
                                                 double deg)
{
    //TODO Documentation
    Image2D::Pointer ret;
    Image2D::SizeType size;
    Image2D::RegionType r;

    size[0] = vm[0]->GetLargestPossibleRegion().GetSize()[0];
    size[1] = vm[0]->GetLargestPossibleRegion().GetSize()[1];

    ret = Image2D::New();
    r.SetSize(size);
    ret->SetRegions(r);
    ret->Allocate();

    ImagePolynomial::Pointer polynomialTable = ConstructPolynomialTable(ip);

    itk::ImageRegionIteratorWithIndex<Image2D> it(vm[0], vm[0]->GetLargestPossibleRegion());
    vnl_vector<double> variable(vm.size());

    for(;!it.IsAtEnd(); ++it)
    {
        for(int i = 0; i < variable.size(); i++)
            variable[i] = (*vm[i])[it.GetIndex()];

        vnl_vector<double> variableCoeff = ConstructRow(variable, deg);
        (*ret)[it.GetIndex()] = (*(polynomialTable))[it.GetIndex()](variableCoeff, true);

    }


    return ret;
}


//Display function
void SpectralFilter::SeeImageValue(Image3D::Pointer p)
{
    Image3D::IndexType ind;
    itk::ImageRegionIteratorWithIndex<Image3D> it(p, p->GetLargestPossibleRegion());
    for(it.GoToBegin(); it.GetIndex()[2] == 0; ++it)
    {
        printf("[%i,%i]", it.GetIndex()[0], it.GetIndex()[1]);
        ind[0] = it.GetIndex()[0];
        ind[1] = it.GetIndex()[1];
        for(int k = 0; k < p->GetLargestPossibleRegion().GetSize()[2]; k++)
        {
            ind[2] = k;

            printf("|%0.50f|",(*p)[ind] );
        }
        printf("\n");
    }

}

void SpectralFilter::ShowEstimatedError(std::vector<Image3D::Pointer> vm, Image3D::Pointer pol,
                                        vnl_vector<double> e, double deg)
{
    Image3D::IndexType index;
    ImagePolynomial::IndexType index2;

    //Initialization
    ImagePolynomial::Pointer polynomialTable = ConstructPolynomialTable(pol);
    itk::ImageRegionIteratorWithIndex<Image3D> it(vm[0], vm[0]->GetLargestPossibleRegion());
    for(; it.GetIndex()[2] == 0; ++it)
    {
            //Index Update
            index[0] = it.GetIndex()[0];
            index[1] = it.GetIndex()[1];
            index2[0] = it.GetIndex()[0];
            index2[1] = it.GetIndex()[1];
            //Construction of k pixel line (size P). with k number of energy level
            std::vector<vnl_vector<double> > vom;
            for(int k = 0; k < vm.size(); k++)
                vom.push_back(TakePixelVector(vm[k], index[0], index[1]));

            //Construction of am matrix of the linear system Ax=b with A = am
            //The am matrix is based on vom values
            vnl_sparse_matrix<double> am = ConstructMatrix(vom, deg);

            std::cout<<"Error "<<index2<<" = "
               <<CalculateEstimatedError(am,(*polynomialTable)[index2].GetCoefficients() , e)<<std::endl;
    }
}

// I/O function
vnl_vector<double> SpectralFilter::OpenThickness(std::string name)
{
    std::vector<double> temp;
    std::ifstream stream;
    stream.open(name.c_str(), std::fstream::in);
    if(!stream.is_open())
    {
        throw std::string("Unable to open " + name);
    }

    double read;

    while(!stream.eof())
    {
        stream>>read;
        temp.push_back(read);
    }

    stream.close();
    vnl_vector<double> ret(temp.size());;

    for(int i = 0; i < temp.size(); i++)
    {
        ret[i] = temp[i];
    }
    std::cout<<"size "<<ret.size()<<std::endl;
    return ret;
}

SpectralFilter::Image3D::Pointer SpectralFilter::OpenImage3D(std::string name)
{
    itk::ImageFileReader<Image3D>::Pointer stream = itk::ImageFileReader<Image3D>::New();
    stream->SetFileName(name);
    stream->Update();
    Image3D::Pointer ret = stream->GetOutput();

    return ret;
}

void SpectralFilter::SaveImage3D(Image3D::Pointer image, std::string name)
{
    itk::ImageFileWriter<Image3D>::Pointer writer = itk::ImageFileWriter<Image3D>::New();
    writer->SetFileName(name);
    writer->SetInput(image);
    writer->Update();
}

SpectralFilter::Image2D::Pointer SpectralFilter::OpenImage2D(std::string name)
{
    itk::ImageFileReader<Image2D>::Pointer stream = itk::ImageFileReader<Image2D>::New();
    stream->SetFileName(name);
    stream->Update();
    Image2D::Pointer ret = stream->GetOutput();

    return ret;
}

void SpectralFilter::SaveImage2D(Image2D::Pointer image, std::string name)
{
    itk::ImageFileWriter<Image2D>::Pointer writer = itk::ImageFileWriter<Image2D>::New();
    writer->SetFileName(name);
    writer->SetInput(image);
    writer->Update();
}

//private function

vnl_vector<double> SpectralFilter::TakePixelVector(Image3D::Pointer image, int i, int j)
{
    Image3D::RegionType region1 = image->GetLargestPossibleRegion();
    vnl_vector<double> ret(region1.GetSize()[2]);
    Image3D::IndexType index;
    index[0] = i;
    index[1] = j;

    //Take a pixel line along the third dimension of the image, pick P pixel
    for(int k = 0; k < region1.GetSize()[2]; k++)
    {
        index[2] = k;
        ret[k] = (*image)[index];
    }

    return ret;
}

vnl_vector<double> SpectralFilter::ConstructRow(vnl_vector<double> variables, double degree)
{
    vnl_vector<double> ret;
    //Return special row for a 2.5 degree polynomial of 2 variables
    if(degree == 2.5 && variables.size() == 2)
    {
        ret.set_size(8);
        ret[0] = 1;
        ret[1] = variables[0]; //m1
        ret[2] = variables[1]; //m2
        ret[3] = variables[0]*variables[0]; //m1²
        ret[4] = variables[0]*variables[1]; //m1*m2<double, 3>
        ret[5] = variables[1]*variables[1]; //m2²
        ret[6] = variables[0]*variables[0]*variables[0]; //m1^3
        ret[7] = variables[1]*variables[1]*variables[1]; //m2^3
    }
    else
        ret = PolynomialCalculator::PolynomialDetermination(variables, degree); //General case

    return ret;
}

vnl_sparse_matrix<double> SpectralFilter::ConstructMatrix(const std::vector<vnl_vector<double> > & vom,
                                                          double deg)
{
    vnl_sparse_matrix<double> ret;
    ret.resize(vom[0].size(), PolynomialSize(vom.size(), deg));

    //Contruction of a Pxl matrix. with l determined by k energy level and degree
    //e.g: For k = 2 and degree = 2. l = 6
    for(int i = 0; i < vom[0].size(); i++)
    {
        vnl_vector<double> vep(vom.size());
        //Construction of variables set (taking i-st value of each vnl_vector)
        for(int j = 0; j < vom.size(); j++)
            vep[j] = vom[j][i];

        vnl_vector<double> temp = ConstructRow(vep,deg);
        //Insert the i-st row of the matrix
        for(int j = 0; j < temp.size(); j++)
            ret.put(i,j,temp[j]);

    }


    return ret;
}

SpectralFilter::ImagePolynomial::Pointer SpectralFilter::ConstructPolynomialTable(Image3D::Pointer ip)
{
    ImagePolynomial::Pointer ret;
    ImagePolynomial::RegionType region;
    ImagePolynomial::SizeType size;
    ImagePolynomial::IndexType index;
    Image3D::IndexType indexP;

    //Initialization of ret
    size[0] = ip->GetLargestPossibleRegion().GetSize()[0];
    size[1] = ip->GetLargestPossibleRegion().GetSize()[1];
    ret = ImagePolynomial::New();
    region.SetSize(size);
    ret->SetRegions(region);
    ret->Allocate();
    itk::ImageRegionIteratorWithIndex<Image3D> it(ip, ip->GetLargestPossibleRegion());
    int p = ip->GetLargestPossibleRegion().GetSize()[2];
    vnl_vector<double> coefficient(p);
    for(; it.GetIndex()[2] == 0; ++it)
    {
        index[0] = it.GetIndex()[0];
        index[1] = it.GetIndex()[1];
        indexP[0] = it.GetIndex()[0];
        indexP[1] = it.GetIndex()[1];

        //Polynomial reading from original matrix
        for(int i = 0; i < p; i++)
        {
            indexP[2] = i;
            coefficient[i] = (*ip)[indexP];
        }
        (*ret)[index].SetCoefficient(coefficient);
    }
    return ret;
}

int SpectralFilter::PolynomialSize(int nbVariables, double degree)
{
    //Return special value for 2.5 degree polynomial of 2 variables
    if(degree == 2.5 && nbVariables == 2)
        return 8;
    else
        return PolynomialCalculator::NumberOfCoefficient(nbVariables, degree); //General case

}

std::vector<double> SpectralFilter::SolveLinearSystem(vnl_sparse_matrix<double> am, vnl_vector<double> ep)
{
    vnl_vector<double> coeficient(am.columns());
    coeficient.fill(0);
    std::vector<double> ret;
    vnl_sparse_matrix_linear_system<double> temp(am, ep);
    vnl_lsqr resolution(temp);
    resolution.set_max_iterations(1000000);
    resolution.minimize(coeficient); 

    for(int i = 0 ; i < am.columns(); i++)
    {
        ret.push_back(coeficient[i]);
    }

    std::cout<<"Error = "<<resolution.get_resid_norm_estimate()<<std::endl<<"Nb iteration"
       <<resolution.get_number_of_iterations()<<std::endl;

    return ret;
}

double SpectralFilter::CalculateEstimatedError(vnl_sparse_matrix<double> a, vnl_vector<double> x,
                                               vnl_vector<double> b)
{
    vnl_vector<double> ax(a.rows());
    a.mult(x, ax);
    vnl_vector<double> axb(ax.size());
    axb = ax-b;

    double ret = 0;

    for(int i = 0; i < axb.size(); i++)
        ret += (axb(i)*axb(i));

    return std::sqrt(ret);
}

} //end namespace rtk

#endif //SPECTRAL_HXX
