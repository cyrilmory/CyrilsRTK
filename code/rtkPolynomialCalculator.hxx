#ifndef POLYNOMIALCALCULATOR_HXX
#define POLYNOMIALCALCULATOR_HXX

#include "rtkPolynomialCalculator.h"

namespace rtk
{
PolynomialCalculator::PolynomialCalculator()
{}

PolynomialCalculator::PolynomialCalculator(vnl_vector<double> coefficients, int degree): m_coefficients(coefficients), m_degree(degree)
{}

vnl_vector<double> PolynomialCalculator::PolynomialDetermination(vnl_vector<double> variables, int degree)
{
    vnl_vector<double> ret (NumberOfCoefficient(variables.size(), degree)); //Return variable


    std::vector<std::vector<int> > start; //Working variable used to performed correct start for each variable at each iteration
    start.push_back(std::vector<int>());

    //Fill the first start with 0
    start[0].resize(variables.size());
    fill(start[0].begin(), start[0].end(), 0);

    std::vector<std::vector<double> > it; //Working variables for stocking results
    it.resize(degree+1);
    it[0].push_back(1);
    int c = 1;
    ret[0] = 1;

    //Compute the i-st degree
    for(int i = 0; i < degree; i ++)
    {
        start.push_back(std::vector<int>());
        start.back().resize(variables.size());

        //Compute the j-st variables during the i-st degree
        for(int j = 0; j < variables.size(); j++)
        {
            //Compute values from the correct start to the last (i-1)-st degree value
            for(int k = start[i][j]; k < it[i].size(); k++, c++)
            {
                it[i+1].push_back(variables[j]*it[i][k]);
                ret[c] = variables[j]*it[i][k];
                if((j+1)<variables.size())
                    start[i+1][j+1]=std::max(it[i+1].size(),(unsigned long)(j+1));
            }
        }
    }


    return ret;

}

double PolynomialCalculator::PolynomialApplication(vnl_vector<double> variables, bool init) const
{
    return (*this)(variables, init);
}

double PolynomialCalculator::operator ()(vnl_vector<double> variablesSet, bool isInit) const
{
    double ret = 0;
    //Compute de indeterminate vector from variablesSet
    vnl_vector<double> variables;

    if(!isInit)
        variables = PolynomialDetermination(variablesSet, GetDegree());
    else
        variables = variablesSet;


    if(variables.size() != GetCoefficients().size())
        throw std::string("Imcompatible vectors size");

    //Compute sum(coefficients[i]*variables[i]
    //coefficients.size() == variables.size();
    for(int i = 0; i < m_coefficients.size(); i++)
        ret += (m_coefficients[i]*variables[i]);

    return ret;
}

int PolynomialCalculator::NumberOfCoefficient(int nbVariable, int degree)
{
    int ret = 0;
    std::vector<int> temp;
    temp.resize(degree+1);

    if(degree == 0)
    {
        ret = 1;
    }
    else
    {
        //Fill the working vector with 1 value
        std::fill(temp.begin(), temp.end(), 1);

        //Compute binomial coefficient
        //Example: for 2 variables the vector compute will be [1,3,6,10,15], i.e 1 coefficient for 0 order polynomial, 3 for 1st order, 6 for 2nd order ....
        for(int i = 1; i < nbVariable+1; i++)
            for(int j = 1; j < temp.size(); j++)
                temp[j] += temp[j-1];
    }

    return temp[degree];
}

vnl_vector<double> PolynomialCalculator::GetCoefficients() const
{
    return m_coefficients;
}

void PolynomialCalculator::SetCoefficient(vnl_vector<double> coefficient)
{
    m_coefficients = coefficient;
}

int PolynomialCalculator::GetDegree() const
{
    return m_degree;
}

void PolynomialCalculator::SetDegree(int degree)
{
    m_degree = degree;
}

}//end namespace rtk

#endif //POLYNOMIALCALCULATOR_HXX
