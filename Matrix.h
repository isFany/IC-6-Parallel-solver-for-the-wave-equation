/**
 * @file Matrix.h
 * @author Fan Yang
 * @version 1.0
 * @date 2021.04.12
 * @brief The header file of Matrix.cpp
 * @details
 * This class is used to create a matrix. This class has six parameters. The matrix rows, columns
 * and values. The other three parameters are used to determine the boundary conditions of the matrix.
 * Some functions that are used to the constructor, print the matrix values, initialize, evolve and build 
 * boundaries are also provided. 
 */


#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;
/**
 * @brief To get the circumference ratio 
 */
#define PI acos(-1)


const int imax = 100;  ///< specify the number of matrix's rows
const int jmax = 100;   ///< specify the number of matrix's rows
const double t_max = 30.0;   ///< specify t_max
// double t = 0;
// double t_out = 0.0;
const double dt_out = 0.1;   ///< The time interval to generate the .dat file
const double y_max = 10.0;  ///<  specify the maximum value of y
const double x_max = 10.0;   ///<  specify the maximum value of x
const double dx = x_max / ((double)imax - 1);  ///<  specify the value of dx
const double dy = y_max / ((double)jmax - 1);   ///<  specify the value of dy
const double c = 1;     ///< specify the value of c in equation. c is the speed of the wave
const double dt = 0.1 * min(dx, dy) / c;    ///< each evolution needs dt

const double r_splash = 1.0;   ///< it is used for initialize
const double x_splash = 3.0;  ///< it is used for initialize
const double y_splash = 3.0;   ///< it is used for initialize


/**
* @class <Matrix>
* @brief This is matrix class to store basic information of matrix
*/
class Matrix {
public:

	double* values = nullptr;  ///< It is a one-dimensional array to store the all matrix value

	int row = -1;  ///< The value of matrix's rows

	int col =-1;  ///< The value of matrix's cols

	bool  periodic = false;  ///< boolean value to specify whether this matrix is periodic or not

	bool Dirichlet = false;  ///< boolean value to specify whether this matrix's boundray is Dirichlet or not

	bool Gradient = false;   ///< boolean value to specify whether this matrix's boundray is Gradient or not
	
	/**
	* @brief		The constructor of Matrix
	* @param[in]	row : The rows of the matrix
	* @param[in]	col : The columns of the matrix
	* @param[in]	periodic : judge whether the boundray is periodic or not
	* @param[in]	Dirichlet : judge whether the boundray is Dirichlet or not
	* @param[in]	Gradient : judge whether the boundray is Neumann or not
	* @details  
	* This is constructor. It is used to create an  object of the matrix. 
	* Because we have a boundary. So add 2 to the number of rows and columns passed in
	* @note	Notice that the value of row and col is bigger than 0.
	* @par Sample
	* @code
	*	Matrix A(10,10,false,true,false);
	*	cout<<A.row<<A.col<<endl;   //Both row and col of A are 12
	* @endcode
	*/
	Matrix(int row, int col, bool  periodic, bool  Dirichlet,bool Gradient);

	/**
	* @brief		print the value of the matrix. It is used to test
	* @return        void
	*/
	void printMatrix();

	/**
	* @brief		This function is used to the boundary of the matrix
	* @return        void
	* @details       
	* When the  boundary is Dirichlet, it means the value of the top and bottom rows, and left
	* and right columns of the matrix are zero. 
	*/
	void built_Dirichlet(int rows,int cols);
	
	/**
    * @brief		This function is used to the boundary of the matrix
    * @return        void
    * @details       
	* When the boundary is Neumann, it means the value of the first row of the matrix 
	* is equal to the value of the second row of the matrix. The value of the last row 
	* of the matrix is equal to the value of the second-to-last row of the matrix. The 
	* value in the first column of the matrix is equal to the value in the second column.
	* The value in the last column of the matrix is equal to the value in the penultimate
	* column.
    */
	void built_Gradient(int rows, int cols);
	
	/**
    * @brief		This function is used to the boundary of the matrix
    * @return       void
    * @details       
	* When the boundary is Periodic, it means the value of the first row  of the matrix is 
	* equal to the value of the last row. The value of the first column of the matrix is equal
	* to the value of the last column
    */
	void built_Periodic(int rows, int cols);

	/**
    * @brief		The initialization of the matrix
	* @param[in]    rows: The rows of the matrix
	* @param[in]    cols: The columns of the matrix  
    * @return       void
    */
	void initial(int rows, int cols);

	/**
    * @brief		The initialization of the matrix
    * @param[in]    new_grid: the pointer of the new_grid whose type is Matrix
    * @param[in]    grid: the pointer of the grid whose type is Matrix
	* @param[in]    old_grid: the pointer of the old_grid whose type is Matrix
	* @param[in]    id: the id of the current process
	* @param[in]    t: The time
    * @return       void
	* @details
	* The evolution is according to the wave equation. After evolution, according to the 
	* different boundray, assign different values to the boundary
    */
	void evolution(Matrix*& new_grid,Matrix*& grid, Matrix*& old_grid,int id, double& t);

	};