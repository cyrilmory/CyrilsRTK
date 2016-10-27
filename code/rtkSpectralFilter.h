#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <cstdlib>
#include <cmath>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "rtkPolynomialCalculator.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_lsqr.h"
#include "vnl/vnl_sparse_matrix_linear_system.h"
#include "vnl/vnl_inverse.h"


namespace rtk
{

    class SpectralFilter
    {
    public:
        typedef itk::Image<double,3> Image3D;
        typedef itk::Image<double,2> Image2D;
        typedef itk::Image<PolynomialCalculator, 2> ImagePolynomial;

        SpectralFilter();

        /**
         * @brief Process spectral decomposition
         * @param thickness
         * @param images: Calibration images, used to create de the matrix A of the linear system Ax=b.
         * @param Degree: polynomials' degree
         * @return Polynomial
         */
        Image3D::Pointer Decomposition(vnl_vector<double> thickness,
              const std::vector<Image3D::Pointer>& images, double degree);

        /**
           * @brief Process polynomial application to images
           * @param polynomials
           * @param images
           * @param degree
           * @return images
           */
          Image2D::Pointer Forward(Image3D::Pointer polynomials, std::vector<Image2D::Pointer> images, double degree);

          /**
         * @brief Display value from a 3D Image
         * @param image
         */
        void SeeImageValue(Image3D::Pointer);

        void ShowEstimatedError(std::vector<Image3D::Pointer>, Image3D::Pointer, vnl_vector<double>, double);


        /**
         * @brief Read thickness values from a tesxt file
         * @param filename
         */
        static vnl_vector<double> OpenThickness(std::string);

        /**
         * @brief Open a 3D image from a file
         * @param filename
         */
        static Image3D::Pointer OpenImage3D(std::string);

        /**
         * @brief Save a 3D image in a file
         * @param Image3D
         * @param filename
         */
        static void SaveImage3D(Image3D::Pointer, std::string);

        /**
         * @brief Open a 2D image from a file
         * @param filename
         */
        static Image2D::Pointer OpenImage2D(std::string);

        /**
         * @brief Save a 2D in a file
         * @param Image2D
         * @param filename
         */
        static void SaveImage2D(Image2D::Pointer, std::string);

    private:
        /**
         * @brief TakePixelVector function return return the pixel's line of coordinate [i,j]
         * @param image
         * @param i
         * @param j
         */
        vnl_vector<double> TakePixelVector(Image3D::Pointer image, int i, int j);

        /**
         * @brief Compute polynomial inderminate. Special case can be implemented instead of general case
         * @param variables: Variables set
         * @param degree
        **/
        vnl_vector<double> ConstructRow(vnl_vector<double> variables, double degree = 2);
        /**
         * @brief Construct the matrix of a linear system Ax=b
         * @param pixelsLines: Vector of line of pixels
         * @param degree
         * @return A
         */
        vnl_sparse_matrix<double> ConstructMatrix(const std::vector<vnl_vector<double> >& pixelsLines, double degree);
        /**
         * @brief Construct a image of PolynomialCalculator from an image packaging polynomials coefficients
         * @param polynomialImage
         */
        ImagePolynomial::Pointer ConstructPolynomialTable(Image3D::Pointer polynomialImage);
        /**
         * @brief PolynomialSize return the number of coefficient from a polynomial of n-degree and k-variables
         * @param nbVariables
         * @param degree
         */
        int PolynomialSize(int nbVariables, double degree);
        /**
         * @brief Solve the linear system Ax=b by a LeastSquare (LSQR) method
         * @param a
         * @param b
         * @return x
         */
        std::vector<double> SolveLinearSystem(vnl_sparse_matrix<double> a, vnl_vector<double> b);
        /**
         * @brief Compute estimated LSQR error.
         * @param a
         * @param x
         * @param b
         * @return ||Ax-b||
         */
        double CalculateEstimatedError(vnl_sparse_matrix<double> a, vnl_vector<double> x, vnl_vector<double> b);
    };
}
#include "rtkSpectralFilter.hxx"

#endif // SPECTRAL_H
