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

#ifndef rtkSchlomka2008NegativeLogLikelihood_h
#define rtkSchlomka2008NegativeLogLikelihood_h

#include "rtkProjectionsDecompositionNegativeLogLikelihood.h"

#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkVariableSizeMatrix.h>

namespace rtk
{
  /** \class rtkSchlomka2008NegativeLogLikelihood
   * \brief Cost function from the Schlomka 2008 PMB paper
   *
   * See the reference paper: "Experimental feasibility of multi-energy photon-counting
   * K-edge imaging in pre-clinical computed tomography", Schlomka et al, PMB 2008
   *
   * \author Cyril Mory
   *
   * \ingroup ReconstructionAlgorithm
   */

// We have to define the cost function first
class Schlomka2008NegativeLogLikelihood : public rtk::ProjectionsDecompositionNegativeLogLikelihood
{
public:

  typedef Schlomka2008NegativeLogLikelihood                     Self;
  typedef rtk::ProjectionsDecompositionNegativeLogLikelihood    Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( Schlomka2008NegativeLogLikelihood, rtk::ProjectionsDecompositionNegativeLogLikelihood );

//  enum { SpaceDimension=m_NumberOfMaterials };

  typedef Superclass::ParametersType          ParametersType;
  typedef Superclass::DerivativeType          DerivativeType;
  typedef Superclass::MeasureType             MeasureType;

  typedef Superclass::DetectorResponseType      DetectorResponseType;
  typedef Superclass::MaterialAttenuationsType  MaterialAttenuationsType;
  typedef Superclass::MeasuredDataType          MeasuredDataType;
  typedef Superclass::IncidentSpectrumType      IncidentSpectrumType;

  // Constructor
  Schlomka2008NegativeLogLikelihood()
  {
  }

  // Destructor
  ~Schlomka2008NegativeLogLikelihood()
  {
  }

  vnl_vector<double> ForwardModel(const ParametersType & lineIntegrals) const ITK_OVERRIDE
  {
  // Variable length vector and variable size matrix cannot be used in linear algebra operations
  // Get their vnl counterparts, which can
  vnl_vector<double> vnl_vec(GetAttenuatedIncidentSpectrum(lineIntegrals).GetDataPointer(), GetAttenuatedIncidentSpectrum(lineIntegrals).GetSize());

  // Apply detector response, getting the lambdas
  return (m_DetectorResponse.GetVnlMatrix() * vnl_vec);
  }

  itk::VariableLengthVector<double> GetAttenuatedIncidentSpectrum(const ParametersType & lineIntegrals) const
  {
  // Solid angle of detector pixel, exposure time and mAs should already be
  // taken into account in the incident spectrum image

  // Apply attenuation at each energy
  itk::VariableLengthVector<float> attenuatedIncidentSpectrum;
  attenuatedIncidentSpectrum.SetSize(m_NumberOfEnergies);
  attenuatedIncidentSpectrum.Fill(0);
  for (unsigned int e=0; e<m_NumberOfEnergies; e++)
    {
    float totalAttenuation = 0.;
    for (unsigned int m=0; m<m_NumberOfMaterials; m++)
      {
      totalAttenuation += lineIntegrals[m] * m_MaterialAttenuations[m][e];
      }

    attenuatedIncidentSpectrum[e] = m_IncidentSpectrum[e] * std::exp(-totalAttenuation);
    }

  return attenuatedIncidentSpectrum;
  }

  itk::VariableLengthVector<float> GetInverseCramerRaoLowerBound(const ParametersType & lineIntegrals) const
  {
  // Get some required data
  vnl_vector<double> attenuatedIncidentSpectrum(GetAttenuatedIncidentSpectrum(lineIntegrals).GetDataPointer(), GetAttenuatedIncidentSpectrum(lineIntegrals).GetSize());
  vnl_vector<double> lambdas = ForwardModel(lineIntegrals);

  // Compute the vector of m_b / lambda_b²
  vnl_vector<double> weights;
  weights.set_size(m_NumberOfSpectralBins);
  for (unsigned int i=0; i<m_NumberOfSpectralBins; i++)
    weights[i] = m_MeasuredData[i] / (lambdas[i] * lambdas[i]);

  // Prepare intermediate variables
  vnl_vector<double> intermediate_a;
  vnl_vector<double> intermediate_a_prime;
  vnl_vector<double> partial_derivative_a;
  vnl_vector<double> partial_derivative_a_prime;

  // Compute the Fischer information matrix
  itk::VariableSizeMatrix<float> Fischer;
  Fischer.SetSize(m_NumberOfMaterials, m_NumberOfMaterials);
  for (unsigned int a=0; a<m_NumberOfMaterials; a++)
    {
    for (unsigned int a_prime=0; a_prime<m_NumberOfMaterials; a_prime++)
      {
      // Compute the partial derivatives of lambda_b with respect to the material line integrals
      intermediate_a = element_product(attenuatedIncidentSpectrum, m_MaterialAttenuations.GetVnlMatrix().get_row(a));
      intermediate_a_prime = element_product(attenuatedIncidentSpectrum, m_MaterialAttenuations.GetVnlMatrix().get_row(a_prime));

      partial_derivative_a = m_DetectorResponse.GetVnlMatrix() * intermediate_a;
      partial_derivative_a_prime = m_DetectorResponse.GetVnlMatrix() * intermediate_a_prime;

      // Multiply them together element-wise, then dot product with the weights
      partial_derivative_a_prime = element_product(partial_derivative_a, partial_derivative_a_prime);
      Fischer[a][a_prime] = dot_product(partial_derivative_a_prime,weights);
      }
    }

  // Invert the Fischer matrix
  itk::VariableLengthVector<double> diag;
  diag.SetSize(m_NumberOfMaterials);
  diag.Fill(0);

  Fischer = Fischer.GetInverse();

  // Return the inverses of the diagonal components (i.e. the inverse variances, to be used directly in WLS reconstruction)
  for (unsigned int mat=0; mat<m_NumberOfMaterials; mat++)
    diag[mat] = 1./Fischer[mat][mat];
  return diag;
  }

  // Not used with a simplex optimizer, but may be useful later
  // for gradient based methods
  void GetDerivative( const ParametersType & lineIntegrals,
                      DerivativeType & derivatives ) const ITK_OVERRIDE
  {
  // Set the size of the derivatives vector
  derivatives.set_size(m_NumberOfMaterials);

  // Get some required data
  vnl_vector<double> attenuatedIncidentSpectrum(GetAttenuatedIncidentSpectrum(lineIntegrals).GetDataPointer(), GetAttenuatedIncidentSpectrum(lineIntegrals).GetSize());
  vnl_vector<double> lambdas = ForwardModel(lineIntegrals);

  // Compute the vector of 1 - m_b / lambda_b
  vnl_vector<double> weights;
  weights.set_size(m_NumberOfSpectralBins);
  for (unsigned int i=0; i<m_NumberOfSpectralBins; i++)
    weights[i] = 1 - (m_MeasuredData[i] / lambdas[i]);

  // Prepare intermediate variables
  vnl_vector<double> intermediate_a;
  vnl_vector<double> partial_derivative_a;

  for (unsigned int a=0; a<m_NumberOfMaterials; a++)
    {
    // Compute the partial derivatives of lambda_b with respect to the material line integrals
    intermediate_a = element_product(-attenuatedIncidentSpectrum, m_MaterialAttenuations.GetVnlMatrix().get_row(a));
    partial_derivative_a = m_DetectorResponse.GetVnlMatrix() * intermediate_a;

    // Multiply them together element-wise, then dot product with the weights
    derivatives[a] = dot_product(partial_derivative_a,weights);
    }
  }

  // Main method
  MeasureType  GetValue( const ParametersType & parameters ) const ITK_OVERRIDE
  {
  // Forward model: compute the expected number of counts in each bin
  vnl_vector<double> lambdas = ForwardModel(parameters);

  // Compute the negative log likelihood from the lambdas
  long double measure = 0;
  for (unsigned int i=0; i<m_NumberOfSpectralBins; i++)
    measure += lambdas[i] - std::log((long double)lambdas[i]) * m_MeasuredData[i];

  return measure;
  }

  itkSetMacro(IncidentSpectrum, IncidentSpectrumType)
  itkGetMacro(IncidentSpectrum, IncidentSpectrumType)

  itkSetMacro(NumberOfSpectralBins, unsigned int)
  itkGetMacro(NumberOfSpectralBins, unsigned int)

protected:
  IncidentSpectrumType        m_IncidentSpectrum;
  unsigned int                m_NumberOfSpectralBins;

private:
  Schlomka2008NegativeLogLikelihood(const Self &); //purposely not implemented
  void operator = (const Self &); //purposely not implemented

};

}// namespace RTK

#endif
