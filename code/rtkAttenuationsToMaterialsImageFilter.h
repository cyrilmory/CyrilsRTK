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

#ifndef rtkAttenuationsToMaterialsImageFilter_h
#define rtkAttenuationsToMaterialsImageFilter_h

#include <itkImageToImageFilter.h>

#include "rtkPolynomialCalculator.h"

#include <vnl_sparse_matrix.h>
#if ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
  #include <itkImageRegionSplitterDirection.h>
#endif

namespace rtk
{
  /** \class AttenuationsToMaterialsImageFilter
   * \brief Decomposes the attenuations into materials,
   * based on a prior polynomial calibration
   *
   * \author Cyril Mory
   *
   * \ingroup ReconstructionAlgorithm
   */

template<typename InputImageType, typename OutputImageType>
class ITK_EXPORT AttenuationsToMaterialsImageFilter :
  public itk::ImageToImageFilter<InputImageType, OutputImageType>
{
public:
  /** Standard class typedefs. */
  typedef AttenuationsToMaterialsImageFilter     Self;
  typedef itk::ImageToImageFilter<InputImageType, OutputImageType>     Superclass;
  typedef itk::SmartPointer<Self>                                      Pointer;
  typedef itk::SmartPointer<const Self>                                ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self)

  /** Runtime information support. */
  itkTypeMacro(AttenuationsToMaterialsImageFilter, ImageToImageFilter)

  /** Set/Get macros for the polynomial's degree */
  itkSetMacro(PolynomialDegree, double)
  itkGetMacro(PolynomialDegree, double)

  /** Set/Get the coefficients of the polynomial image */
  void SetInputPolynomialCoefficients(const InputImageType* PolynomialCoefficients);
  typename InputImageType::ConstPointer GetInputPolynomialCoefficients();

protected:
  AttenuationsToMaterialsImageFilter();
  ~AttenuationsToMaterialsImageFilter() ITK_OVERRIDE {}

  void UpdateNumberOfCoefficients();

  void GenerateOutputInformation() ITK_OVERRIDE;

  void GenerateInputRequestedRegion() ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void ThreadedGenerateData(const typename OutputImageType::RegionType& outputRegionForThread, itk::ThreadIdType itkNotUsed(threadId)) ITK_OVERRIDE;

#if ITK_VERSION_MAJOR > 4 || (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 4)
    /** Splits the OutputRequestedRegion along the first direction, not the last */
    const itk::ImageRegionSplitterBase* GetImageRegionSplitter(void) const ITK_OVERRIDE;
    itk::ImageRegionSplitterDirection::Pointer  m_Splitter;
#endif

  /** The inputs should not be in the same space so there is nothing
   * to verify. */
  void VerifyInputInformation() ITK_OVERRIDE {}

  vnl_vector<double> ConstructRow(vnl_vector<double> variables);

  /** Degree of the polynomial. Degree 2.5 corresponds to a special
   * case of degree 3, with some coefficients set to zero.
   * Set by the user. */
  double m_PolynomialDegree;

  /** Number of variables of the polynomial, equal to the number
   * of energy spectra or detector responses, i.e. 2 in dual energy,
   * but the implementation can handle more (e.g. for spectral CT).
   * Determined by the vector length of the input image */
  unsigned int m_NumberOfVariables;

  /** Number of materials used in the calibration process.
   * Often equal to 2 in dual energy, but can be more.
   * Determined from the vector length of the polymial coefficients image */
  unsigned int m_NumberOfMaterials;

  /** Number of coefficients of the polynomial.
   * Computed from the polynomial degree and the number of variables */
  unsigned int m_NumberOfCoefficients;

private:
  //purposely not implemented
  AttenuationsToMaterialsImageFilter(const Self&);
  void operator=(const Self&);

}; // end of class

} // end namespace rtk


#ifndef ITK_MANUAL_INSTANTIATION
#include "rtkAttenuationsToMaterialsImageFilter.hxx"
#endif

#endif
