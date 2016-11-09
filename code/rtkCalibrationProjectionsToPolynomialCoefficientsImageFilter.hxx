/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef rtkCalibrationProjectionsToPolynomialCoefficientsImageFilter_hxx
#define rtkCalibrationProjectionsToPolynomialCoefficientsImageFilter_hxx

#include "rtkCalibrationProjectionsToPolynomialCoefficientsImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableLengthVector.h>

#include <vnl_sparse_matrix_linear_system.h>
#include <vnl_lsqr.h>

namespace rtk
{

template<typename InputImageType, typename OutputImageType>
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::CalibrationProjectionsToPolynomialCoefficientsImageFilter()
{
  m_PolynomialDegree = 2;
  m_NumberOfVariables = 2;
  m_NumberOfMaterials = 2;
  m_NumberOfCoefficients = 6;
  m_NumberOfThicknessesCombinations = 100;

#if ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
  // Set the direction along which the output requested region should NOT be split
  m_Splitter = itk::ImageRegionSplitterDirection::New();
  m_Splitter->SetDirection(OutputImageType::ImageDimension - 1);
#else
  // Old versions of ITK (before 4.4) do not have the ImageRegionSplitterDirection
  // and should run this filter with only one thread
  this->SetNumberOfThreads(1);
#endif
}

template<typename InputImageType, typename OutputImageType>
void
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::UpdateNumberOfCoefficients()
{
  //Return special value for 2.5 degree polynomial of 2 variables
  if(m_PolynomialDegree == 2.5 && m_NumberOfVariables == 2)
    m_NumberOfCoefficients = 8;
  else
    m_NumberOfCoefficients = PolynomialCalculator::NumberOfCoefficient(m_NumberOfVariables, m_PolynomialDegree); //General case
}

template<typename InputImageType, typename OutputImageType>
void
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::GenerateOutputInformation()
{
  m_NumberOfVariables = this->GetInput()->GetVectorLength(); //Should be two for dual energy
  this->UpdateNumberOfCoefficients();
  m_NumberOfThicknessesCombinations = this->GetInput()->GetLargestPossibleRegion().GetSize(InputImageType::ImageDimension - 1);

  typename OutputImageType::RegionType largest = this->GetInput()->GetLargestPossibleRegion();
  largest.SetSize(OutputImageType::ImageDimension - 1, m_NumberOfCoefficients);
  largest.SetIndex(OutputImageType::ImageDimension - 1, 0);

  this->GetOutput()->SetLargestPossibleRegion(largest);
  this->GetOutput()->SetVectorLength(m_NumberOfMaterials);
}

template<typename InputImageType, typename OutputImageType>
void
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::GenerateInputRequestedRegion()
{
  // Input image pointer
  typename InputImageType::Pointer inputPtr =
    const_cast< InputImageType * >( this->GetInput() );
  if ( !inputPtr ) return;

  typename InputImageType::RegionType requested = this->GetOutput()->GetRequestedRegion();
  requested.SetSize(InputImageType::ImageDimension - 1, this->GetInput()->GetLargestPossibleRegion().GetSize(InputImageType::ImageDimension - 1));

  inputPtr->SetRequestedRegion( requested );
}

template<typename InputImageType, typename OutputImageType>
vnl_vector<double>
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::ConstructRow(vnl_vector<double> variables)
{
  vnl_vector<double> systemMatrixRow;
  //Return special row for a 2.5 degree polynomial of 2 variables
  if(m_PolynomialDegree == 2.5 && variables.size() == 2)
    {
    systemMatrixRow.set_size(8);
    systemMatrixRow[0] = 1;
    systemMatrixRow[1] = variables[0]; //m1
    systemMatrixRow[2] = variables[1]; //m2
    systemMatrixRow[3] = variables[0]*variables[0]; //m1²
    systemMatrixRow[4] = variables[0]*variables[1]; //m1*m2<double, 3>
    systemMatrixRow[5] = variables[1]*variables[1]; //m2²
    systemMatrixRow[6] = variables[0]*variables[0]*variables[0]; //m1^3
    systemMatrixRow[7] = variables[1]*variables[1]*variables[1]; //m2^3
    }
  else
    systemMatrixRow = PolynomialCalculator::PolynomialDetermination(variables, m_PolynomialDegree); //General case

  return systemMatrixRow;
}

template<typename InputImageType, typename OutputImageType>
vnl_sparse_matrix<double>
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::ConstructMatrix(const std::vector<vnl_vector<double> > measuredData)
{
  vnl_sparse_matrix<double> systemMatrix;
  systemMatrix.resize(m_NumberOfThicknessesCombinations, m_NumberOfCoefficients);

  //Contruction of a Pxl matrix. with l determined by k energy level and degree
  //e.g: For k = 2 and degree = 2. l = 6
  for(int row = 0; row < m_NumberOfThicknessesCombinations; row++)
    {
    vnl_vector<double> temp = ConstructRow(measuredData[row]);
    //Insert the i-st row of the matrix
    for(int col = 0; col < temp.size(); col++)
        systemMatrix.put(row, col, temp[col]);
    }

    return systemMatrix;
}

#if ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
template<typename InputImageType, typename OutputImageType>
const itk::ImageRegionSplitterBase*
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::GetImageRegionSplitter(void) const
{
  return m_Splitter;
}
#endif

template<typename InputImageType, typename OutputImageType>
void
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::BeforeThreadedGenerateData()
{
#if !(ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4))
  if (this->GetNumberOfThreads() > 1)
    {
    itkWarningMacro(<< "This filter cannot use multiple threads with ITK versions older than v4.4. Reverting to single thread behavior");
    this->SetNumberOfThreads(1);
    }
#endif
}

template<typename InputImageType, typename OutputImageType>
void
CalibrationProjectionsToPolynomialCoefficientsImageFilter<InputImageType, OutputImageType>
::ThreadedGenerateData(const typename OutputImageType::RegionType& outputRegionForThread, itk::ThreadIdType itkNotUsed(threadId))
{
  typename InputImageType::RegionType firstSlice = outputRegionForThread;
  firstSlice.SetSize(InputImageType::ImageDimension - 1, 1);

  // Walk the output's first slice. For each pixel, solve a linear system to determine the coefficients of the polynomials
  itk::ImageRegionIteratorWithIndex<OutputImageType> outputIt (this->GetOutput(), firstSlice);
  itk::ImageRegionConstIteratorWithIndex<InputImageType> inputIt (this->GetInput(), firstSlice);

  while(!outputIt.IsAtEnd())
    {
    // Create an std::vector of vnl_vectors, to store the t-uples of calibration measurements
    std::vector< vnl_vector<double> > measuredData;

    // Read the calibration measurements, store them into the vector
    for (int i=0; i < m_NumberOfThicknessesCombinations; i++)
      {
      typename InputImageType::IndexType inputIndex = inputIt.GetIndex();
      inputIndex[InputImageType::ImageDimension - 1] = i;
      itk::VariableLengthVector<double> inputPixel = this->GetInput()->GetPixel(inputIndex);
      vnl_vector<double> tuple(inputPixel.GetDataPointer(),inputPixel.GetSize());
      measuredData.push_back(tuple);
      }

    // Construct the system matrix
    vnl_sparse_matrix<double> systemMatrix = ConstructMatrix(measuredData);

    // Solve the linear system. The class vnl_sparse_matrix_linear_system can only solve systems
    // of the for Ax=b, with b a vector (not a matrix), so we need to run as many resolutions as
    // the number of materials. For each material, we get a vector of coefficients, which we
    // store into a matrix
    vnl_vector<double> coefficientsVector(m_NumberOfCoefficients);
    vnl_matrix<double> coefficientsMatrix(m_NumberOfCoefficients, m_NumberOfMaterials);

    for (unsigned int col=0; col<m_NumberOfMaterials; col++)
      {
      // Initialize the vector
      coefficientsVector.fill(0);

      // Solve the system and get the result in coefficientsVector
      // We have to create an intermediate vnl_vector from the vnl_matrix's column
      // and then pass it to oneMaterialLinearSystem, otherwise the first two elements
      // are set to almost zero. This looks like an ITK bug.
      vnl_vector<double> matrixColumn = m_ThicknessesMatrix.get_column(col);
      vnl_sparse_matrix_linear_system<double> oneMaterialLinearSystem(systemMatrix, matrixColumn);
      vnl_lsqr solver(oneMaterialLinearSystem);
      solver.set_max_iterations(1000000);
      solver.minimize(coefficientsVector);

      // Display information about the convergence
//      std::cout << "Error = " << solver.get_resid_norm_estimate() << std::endl
//                << "Nb iterations = " << solver.get_number_of_iterations() << std::endl;

      // Store coefficientsVector in coefficientsMatrix
      coefficientsMatrix.set_column(col, coefficientsVector);
      }

    // Write into the output image
    itk::VariableLengthVector<double> outputPixel;
    outputPixel.SetSize(m_NumberOfMaterials);
    for (int i=0; i < m_NumberOfCoefficients; i++)
      {
      typename OutputImageType::IndexType outputIndex = outputIt.GetIndex();
      outputIndex[OutputImageType::ImageDimension - 1] = i;

      // Getting coefficientsMatrix.get_row(i).data_block directly results in the
      // first element of the vector being set to zero. This looks like an ITK bug.
      // We have to create an intermediate vnl_vector from the vnl_matrix's row
      // and then get its data_block.
      vnl_vector<double> matrixRow = coefficientsMatrix.get_row(i);
      outputPixel.SetData(matrixRow.data_block(), m_NumberOfMaterials, false);
      this->GetOutput()->SetPixel(outputIndex, outputPixel);
      }

    // Move to next pixel
    ++inputIt;
    ++outputIt;
    }
}

} // end namespace rtk

#endif // rtkCalibrationProjectionsToPolynomialCoefficientsImageFilter_hxx
