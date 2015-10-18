#include "Matrix.h"

Matrix::Matrix(int numRows, int numCols) {
	numRow=numRows;
	numCol=numCols;
	primitive=std::vector< std::vector<double> >(numRow, vec(numCol, 0));
}
Matrix::Matrix(mat &primit) {
	primitive=primit;
	numRow=primit.size();
	numCol=primit[0].size();
}
Matrix::Matrix(vec &primit) { //very inefficient, try to use another constructor if possible
	numRow=primit.size();
	numCol=1;
	primitive=std::vector< std::vector<double> >(numRow, vec(numCol));
	for(int i=0; i<numRow; i++){
		primitive[i][0]=primit[i];
	}
}
std::string Matrix::type(){ //introspection
    return "Matrix";
}
Matrix Matrix::multiply(const Matrix &c) const {
	if(numCol==c.getNumRow()){
		mat retMat(numRow, vec(c.getNumCol(), 0));
		for(int k=0; k<numRow; k++){
			for(int j=0; j<c.getNumCol(); j++){
				for(int i=0; i<numCol; i++){
					retMat[k][j]+=primitive[k][i]*c[i][j];
				}
			}
		}
		return Matrix(retMat);
	}
	else {
		throw std::invalid_argument( "matrix does not match" );
	}
}
Matrix Matrix::multiply(double c) const {
	//Matrix plac=Matrix(c*real, c*im);
	mat retMat(numRow, vec(numCol, 0));
	for(int i=0; i<numRow; i++){
		for(int j=0; j<numCol; j++){
			retMat[i][j]=primitive[i][j]*c;
		}
	}
	return Matrix(retMat);
}
Matrix Matrix::multiply(const vec &c) const {
	if(numCol==c.size()){
		mat retMat(numRow, vec(1, 0));
		for(int k=0; k<numRow; k++){
			for(int i=0; i<numCol; i++){
				retMat[k][0]+=primitive[k][i]*c[i];
			}

		}
		return Matrix(retMat);
	}
	else if(numRow==c.size()){
		mat retMat(1, vec(numCol, 0));
		for(int k=0; k<numRow; k++){
			for(int i=0; i<numCol; i++){
				retMat[0][i]+=primitive[k][i]*c[k];
			}
		}
		return Matrix(retMat);
	}
	else {
		throw std::invalid_argument( "matrix does not match" );
	}
}
Matrix Matrix::add(const Matrix &c) const {
	if(numRow==c.getNumRow()&&numCol==c.getNumCol()){
		mat retMat(numRow, vec(numCol, 0));
		for(int i=0; i<numRow; i++){
			for(int j=0; j<numCol; j++){
				retMat[i][j]=primitive[i][j]+c[i][j];
			}
		}
	}
	else {
		throw std::invalid_argument( "matrix does not match" );
	}
}
Matrix Matrix::add(double c) const {
	mat retMat(numRow, vec(numCol, 0));
	for(int i=0; i<numRow; i++){
		for(int j=0; j<numCol; j++){
			retMat[i][j]=primitive[i][j]+c;
		}
	}
	return Matrix(retMat);
}
Matrix Matrix::subtract(const Matrix &c) const {
	if(numRow==c.getNumRow()&&numCol==c.getNumCol()){
		mat retMat(numRow, vec(numCol, 0));
		for(int i=0; i<numRow; i++){
			for(int j=0; j<numCol; j++){
				retMat[i][j]=primitive[i][j]-c[i][j];
			}
		}
	}
	else {
		throw std::invalid_argument( "matrix does not match" );
	}
}
Matrix Matrix::subtract(double c) const {
	mat retMat(numRow, std::vector<double>(numCol, 0));
	for(int i=0; i<numRow; i++){
		for(int j=0; j<numCol; j++){
			retMat[i][j]=primitive[i][j]-c;
		}
	}
	return Matrix(retMat);
}
int Matrix::getNumCol() const {
	return numCol;
}
double Matrix::l1norm() {
	double den=(double)numCol*numRow;
	double l1norm=0;
	for(int i=0; i<numRow; i++){
		for(int j=0; j<numCol; j++){
			l1norm+=primitive[i][j];
		}
	}
	return l1norm/den;
}
double Matrix::l2norm() {
	double den=(double)numCol*numRow;
	double l2norm=0;
	for(int i=0; i<numRow; i++){
		for(int j=0; j<numCol; j++){
			l2norm+=primitive[i][j]*primitive[i][j];
		}
	}
	return l2norm/den;
}
int Matrix::size() const {
	return numRow;
}
int Matrix::getNumRow() const{
	return numRow;
}
vec Matrix::operator[](const int index) const {
		assert(index >= 0 && index < numRow);
    return primitive[index];
}

Matrix operator+(const Matrix &c1, const Matrix &c2)
{
    return c1.add(c2);
}
Matrix operator+(double c1, const Matrix &c2)
{
    return c2.add(c1);
}
Matrix operator+(const Matrix &c1, double c2)
{
    return c1.add(c2);
}

Matrix operator-(const Matrix &c1, const Matrix &c2)
{
    return c1.subtract(c2);
}
Matrix operator-(const Matrix &c1, double c2)
{
    return c1.subtract(c2);
}
Matrix operator*(const Matrix &c1, const Matrix &c2)
{
    return c1.multiply(c2);
}
Matrix operator*(double c1, const Matrix &c2)
{
    return c2.multiply(c1);
}
Matrix operator*(const Matrix &c1, double c2)
{
    return c1.multiply(c2);
}
Matrix operator*(const Matrix &c1, const vec &c2)
{
    return c1.multiply(c2);
}
Matrix operator*(const vec &c2, const Matrix &c1)
{
    return c1.multiply(c2);
}
