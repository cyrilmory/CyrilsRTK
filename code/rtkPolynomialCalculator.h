#ifndef POLYNOMIALCALCULATOR_H
#define POLYNOMIALCALCULATOR_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "vnl/vnl_vector.h"

namespace rtk
{
    class PolynomialCalculator
    {
    public:
        PolynomialCalculator();

        /**
         * @brief Constructor
         * @param coefficient: Constant coefficient of a polynomial
         * @param degree: degree of a polynomial
         */
        PolynomialCalculator(vnl_vector<double> coefficient, int degree = 2);

        /**
         * @brief Compute polynomial inderminate for a n-degrees polynomial of k-variables
         * @param variables: set of k-variables
         * @param degree
         * @return
         */
        static vnl_vector<double> PolynomialDetermination(vnl_vector<double> variables, int degree = 2);
        /**
        /** * @brief Same as operator()
        /**
        /** */
        double PolynomialApplication(vnl_vector<double>, bool = false) const;

        /**
         * @brief Calculate result of a polynomial from variables
         *  Example: 2nd degree and 2 variables polynomial
         *           PolynomialCalculator polynomial;
         *           vnl_vector<double> coefficient(6);
         *           //fill coefficient
         *           vnl_vector<double> variables(2);
         *           //fill variables
         *           polynomial.setCoefficient(coefficients);
         *           polynomial.setDegree(2);
         *
         *           std::cout<<"Result = "<<p(variables)<<std::endl;
         *
         * @param variables
         *
         */
        double operator()(vnl_vector<double> variables, bool = false) const;

        /**
         * @brief Compute the number of coefficient for a n-degrees polynomial of k-variables
         *  Example: int degree; int nbVariables;
         *           PolynomialCalculator polynomial;
         *           std::cout<<"The number of coefficient for a <<degree<<" polynomal and "<<
         *            nbVariables<<" is "<<polynomial.NumberOfCoefficient(nbVariables, degree)<<std::endl;
         * @param numberOfVariables
         * @param degree
         */
        static int NumberOfCoefficient(int numberOfVariables, int degree);

        /*Getters and setters*/
        vnl_vector<double> GetCoefficients() const;
        void SetCoefficient(vnl_vector<double>);
        int GetDegree() const;
        void SetDegree(int);

    private:
        vnl_vector<double> m_coefficients;
        int m_degree;
    };
 }
#include "rtkPolynomialCalculator.hxx" 
 
#endif // POLYNOMIALCALCULATOR_H
