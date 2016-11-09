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

#ifndef rtkAttenuationsToMaterialsImageFilter_hxx
#define rtkAttenuationsToMaterialsImageFilter_hxx

#include "rtkAttenuationsToMaterialsImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableLengthVector.h>

namespace rtk
{

template<typename InputImageType, typename OutputImageType>
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::AttenuationsToMaterialsImageFilter()
{
  m_PolynomialDegree = 2;
  m_NumberOfVariables = 2;
  m_NumberOfMaterials = 2;
  m_NumberOfCoefficients = 6;

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
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::SetInputPolynomialCoefficients(const InputImageType* PolynomialCoefficients)
{
  this->SetInput("PolynomialCoefficients", const_cast<InputImageType*>(PolynomialCoefficients));
}

template<typename InputImageType, typename OutputImageType>
typename InputImageType::ConstPointer
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::GetInputPolynomialCoefficients()
{
  return static_cast< const InputImageType * >
          ( this->itk::ProcessObject::GetInput("PolynomialCoefficients") );
}

template<typename InputImageType, typename OutputImageType>
void
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
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
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::GenerateOutputInformation()
{
  m_NumberOfVariables = this->GetInput(0)->GetVectorLength(); //Should be two for dual energy
  m_NumberOfMaterials = this->GetInputPolynomialCoefficients()->GetVectorLength(); //Should be two for dual energy
  this->UpdateNumberOfCoefficients();

  if (m_NumberOfCoefficients != this->GetInputPolynomialCoefficients()->GetLargestPossibleRegion().GetSize(InputImageType::ImageDimension - 1))
    itkGenericExceptionMacro( << "Number of polynomial coefficient computed from m_NumberOfVariables and m_PolynomialDegree (i.e. "
                              << m_NumberOfCoefficients << ") is inconsistent with input polynomial coefficients image size on last dimension (i.e. "
                              << this->GetInputPolynomialCoefficients()->GetLargestPossibleRegion().GetSize(InputImageType::ImageDimension - 1) )

  typename OutputImageType::RegionType largest = this->GetInput()->GetLargestPossibleRegion();

  this->GetOutput()->SetLargestPossibleRegion(largest);
  this->GetOutput()->SetVectorLength(m_NumberOfMaterials);
}

template<typename InputImageType, typename OutputImageType>
void
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::GenerateInputRequestedRegion()
{
  // Input projections
  typename InputImageType::Pointer inputPtr0 =
    const_cast< InputImageType * >( this->GetInput(0) );
  if ( !inputPtr0 ) return;

  typename InputImageType::RegionType requested = this->GetOutput()->GetRequestedRegion();
  inputPtr0->SetRequestedRegion( requested );

  // Polynomial coefficients
  typename InputImageType::Pointer inputPtrCoeffs =
    const_cast< InputImageType * >( this->GetInputPolynomialCoefficients().GetPointer() );
  if ( !inputPtrCoeffs ) return;

  requested.SetSize(InputImageType::ImageDimension -1, inputPtrCoeffs->GetLargestPossibleRegion().GetSize(InputImageType::ImageDimension -1));
  inputPtrCoeffs->SetRequestedRegion( requested );
}

template<typename InputImageType, typename OutputImageType>
vnl_vector<double>
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
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

#if ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
template<typename InputImageType, typename OutputImageType>
const itk::ImageRegionSplitterBase*
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::GetImageRegionSplitter(void) const
{
  return m_Splitter;
}
#endif

template<typename InputImageType, typename OutputImageType>
void
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
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
AttenuationsToMaterialsImageFilter<InputImageType, OutputImageType>
::ThreadedGenerateData(const typename OutputImageType::RegionType& outputRegionForThread, itk::ThreadIdType itkNotUsed(threadId))
{


  // Walk the output. For each pixel, evaluate the precomputed polynomials on the attenuations
  // measured for each energy, which estimates the thickness of each material traversed
  itk::ImageRegionIterator<OutputImageType> outputIt (this->GetOutput(), outputRegionForThread);
  itk::ImageRegionConstIteratorWithIndex<InputImageType> inputIt (this->GetInput(0), outputRegionForThread);

  while(!outputIt.IsAtEnd())
    {
    // Read the measured attenuations for the current pixel
    itk::VariableLengthVector<double> attenuations = inputIt.Get();
    vnl_vector<double> vnl_attenuations(attenuations.GetDataPointer(),attenuations.GetSize());

    // Compute the (not yet coefficiented) monomials for this pixel's measured attenuations
    vnl_vector<double> monomials = ConstructRow(vnl_attenuations);

    // Read the polynomials' coefficients
    vnl_matrix<double> coefficientsMatrix(m_NumberOfMaterials, m_NumberOfCoefficients);
    for (unsigned int i=0; i < m_NumberOfCoefficients; i++)
      {
      typename InputImageType::IndexType coeffsIndex = inputIt.GetIndex();
      coeffsIndex[InputImageType::ImageDimension - 1] = i;
      itk::VariableLengthVector<double> coeffsPixel = this->GetInputPolynomialCoefficients()->GetPixel(coeffsIndex);
      vnl_vector<double> tuple(coeffsPixel.GetDataPointer(),coeffsPixel.GetSize());
      coefficientsMatrix.set_column(i, tuple);
      }

    // Multiply the monomials by their respective coefficients,
    // and sum them up (which evaluates the polynomial)
    vnl_vector<double> evaluations = coefficientsMatrix * monomials;

    // Write into the output image
    itk::VariableLengthVector<double> outputPixel;
    outputPixel.SetSize(m_NumberOfMaterials);
    outputPixel.SetData(evaluations.data_block(), m_NumberOfMaterials, false);
    outputIt.Set(outputPixel);

    // Move to next pixel
    ++inputIt;
    ++outputIt;
    }
}

} // end namespace rtk

#endif // rtkAttenuationsToMaterialsImageFilter_hxx
