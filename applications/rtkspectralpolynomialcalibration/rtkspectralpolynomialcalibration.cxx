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

#include "rtkspectralpolynomialcalibration_ggo.h"
#include "rtkGgoFunctions.h"
#include "rtkConfiguration.h"
#include "rtkMacro.h"
#include "rtkGeneralPurposeFunctions.h"
#include "rtkCalibrationProjectionsToPolynomialCoefficientsImageFilter.h"
#include "rtkImageToVectorImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(rtkspectralpolynomialcalibration, args_info);

  typedef float PixelValueType;
  const unsigned int Dimension = 3;

  typedef itk::Image< PixelValueType, Dimension > CalibrationProjectionType;
  typedef itk::VectorImage< PixelValueType, Dimension > VectorCalibrationProjectionType;
  typedef itk::ImageFileReader<CalibrationProjectionType> CalibrationProjectionReaderType;

  typedef rtk::ImageToVectorImageFilter<CalibrationProjectionType, VectorCalibrationProjectionType> ImageToVectorImageType;

  typedef VectorCalibrationProjectionType PolynomialCoefficientsImageType;
  typedef itk::ImageFileWriter<PolynomialCoefficientsImageType> PolynomialCoefficientWriterType;

  // Read material thicknesses
  std::vector<std::vector<double> > thicknessesTable = rtk::ReadDoubleTableFile(args_info.thicknesses_arg);
  vnl_matrix<double> thicknessesMatrix (thicknessesTable.size(), thicknessesTable[0].size());
  for (unsigned int row=0; row < thicknessesMatrix.rows(); row++)
    {
    for (unsigned int col=0; col < thicknessesMatrix.cols(); col++)
      {
      thicknessesMatrix.put(row, col, thicknessesTable[row][col]);
      }
    }

  // Read input calibration projections
  CalibrationProjectionReaderType::Pointer calibrationProjectionReader = CalibrationProjectionReaderType::New();
  calibrationProjectionReader->SetFileName( args_info.input_arg );
  calibrationProjectionReader->Update();

  // Convert the input image into a vector image
  ImageToVectorImageType::Pointer imageToVectorImage = ImageToVectorImageType::New();
  imageToVectorImage->SetInput(calibrationProjectionReader->GetOutput());
  unsigned int NumberOfEnergies = calibrationProjectionReader->GetOutput()->GetLargestPossibleRegion().GetSize()[Dimension - 1] / thicknessesMatrix.rows();
  imageToVectorImage->SetNumberOfChannels(NumberOfEnergies);

  // Create and set the filter
  typedef rtk::CalibrationProjectionsToPolynomialCoefficientsImageFilter<VectorCalibrationProjectionType, PolynomialCoefficientsImageType> CalibrationFilterType;
  CalibrationFilterType::Pointer calibration = CalibrationFilterType::New();
  calibration->SetInput(imageToVectorImage->GetOutput());
  calibration->SetPolynomialDegree(args_info.degree_arg);
  calibration->SetThicknessesMatrix(thicknessesMatrix);

  TRY_AND_EXIT_ON_ITK_EXCEPTION(calibration->Update())

  // Write output
  PolynomialCoefficientWriterType::Pointer writer = PolynomialCoefficientWriterType::New();
  writer->SetInput(calibration->GetOutput());
  writer->SetFileName(args_info.output_arg);
  writer->Update();

  return EXIT_SUCCESS;
}
