#ifndef __MATRIX_H_INCLUDED__
#define __MATRIX_H_INCLUDED__

#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cassert>
//#include "Number.h"
typedef std::vector<double> vec;
typedef std::vector<vec> mat;
class Matrix {
	private:
		int numRow;
		int numCol;
		mat primitive;
	public:

		friend Matrix operator+(const Matrix &c1, const Matrix &c2) ;
		friend Matrix operator+(const double c1, const Matrix &c2) ;
		friend Matrix operator+(const Matrix &c1, double c2) ;
		friend Matrix operator-(const Matrix &c1, const Matrix &c2) ;
		//friend Matrix operator-(const double c1, const Matrix &c2) ;
		friend Matrix operator-(const Matrix &c1, double c2) ;
		friend Matrix operator*(const Matrix &c1, const Matrix &c2) ;
		friend Matrix operator*(double c1, const Matrix &c2);
		friend Matrix operator*(const Matrix &c1, double c2);
		friend Matrix operator*(const Matrix &c1, const vec &c2); //
		friend Matrix operator*(const vec &c1, const Matrix &c2);//

	//	Matrix(vector<double>&, bool);
		Matrix(mat&);
		Matrix(vec&);
		Matrix(int, int); //empty
		Matrix multiply(const Matrix&) const;
		double l1norm();
		double l2norm();
		std::string type();
		Matrix multiply(double) const;
		Matrix multiply(const std::vector<double>&) const;
		Matrix add(const Matrix&) const;
		Matrix add(double) const;
		Matrix subtract(const Matrix&) const;
		Matrix subtract(double) const;
		vec operator[](const int index) const;//
	//	Matrix power(double) const;
	//Matrix divide(const Complex&) const;
		//Matrix divide(double) const;

		//double getPrimitive() const;
		int getNumCol() const;
		int getNumRow() const;
		int size() const;
		//double getIm() const;

};



#endif
