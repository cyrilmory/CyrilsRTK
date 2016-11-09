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

#include "rtkattenuationstomaterials_ggo.h"
#include "rtkGgoFunctions.h"
#include "rtkConfiguration.h"
#include "rtkMacro.h"
#include "rtkGeneralPurposeFunctions.h"
#include "rtkAttenuationsToMaterialsImageFilter.h"
#include "rtkImageToVectorImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(rtkattenuationstomaterials, args_info);

  typedef float PixelValueType;
  const unsigned int Dimension = 3;

  typedef itk::Image< PixelValueType, Dimension > ProjectionStackType;
  typedef itk::VectorImage< PixelValueType, Dimension > VectorProjectionStackType;

  typedef itk::ImageFileReader<ProjectionStackType> ProjectionStackReaderType;
  typedef itk::ImageFileReader<VectorProjectionStackType> VectorProjectionStackReaderType;
  typedef rtk::ImageToVectorImageFilter<ProjectionStackType, VectorProjectionStackType> ImageToVectorImageType;
  typedef itk::ImageFileWriter<VectorProjectionStackType> MaterialDecomposedProjectionStackWriterType;

  // Read input measured projections
  ProjectionStackReaderType::Pointer projectionsReader = ProjectionStackReaderType::New();
  projectionsReader->SetFileName( args_info.input_arg );
  projectionsReader->Update();

  // Convert the input image into a vector image
  ImageToVectorImageType::Pointer imageToVectorImage = ImageToVectorImageType::New();
  imageToVectorImage->SetInput(projectionsReader->GetOutput());
  imageToVectorImage->SetNumberOfChannels(args_info.nenergies_arg);

  // Read input measured projections
  VectorProjectionStackReaderType::Pointer polynomialsCoefficientsReader = VectorProjectionStackReaderType::New();
  polynomialsCoefficientsReader->SetFileName( args_info.coeffs_arg );
  polynomialsCoefficientsReader->Update();

  // Create and set the filter
  typedef rtk::AttenuationsToMaterialsImageFilter<VectorProjectionStackType, VectorProjectionStackType> AttenuationsToMaterialsType;
  AttenuationsToMaterialsType::Pointer forward = AttenuationsToMaterialsType::New();
  forward->SetInput(0, imageToVectorImage->GetOutput());
  forward->SetInputPolynomialCoefficients(polynomialsCoefficientsReader->GetOutput());
  forward->SetPolynomialDegree(args_info.degree_arg);

  TRY_AND_EXIT_ON_ITK_EXCEPTION(forward->Update())

  // Write output
  MaterialDecomposedProjectionStackWriterType::Pointer writer = MaterialDecomposedProjectionStackWriterType::New();
  writer->SetInput(forward->GetOutput());
  writer->SetFileName(args_info.output_arg);
  writer->Update();

  return EXIT_SUCCESS;
}
