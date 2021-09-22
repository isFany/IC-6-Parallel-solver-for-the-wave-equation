/**
 * @file Solver.h
 * @author Fan Yang
 * @version 1.0
 * @date 2021.04.12
 * @brief The header file of Solver.cpp
 * @details
 * This file provides two Matrix subclasses called strip and grid. These two classes solve 
 * the wave equation using two decompositions, stripe decomposition and grid decomposition. 
 * These functions are using MPI to do the parallel. By creating dateType, only the marginal
 * value of the matrix is sent. In addition, also provide the function to write the matrix value
 * to .dat file. 
 */

#pragma once

#include "Matrix.h"
#include <mpi.h>
#include <algorithm> 

 /**
 * @class <strip>
 * @brief 
 * The class strip is the subclass of the Matrix. The object of this class can solve wave function using 
 * stripe decomposition.
 */
class strip : public Matrix
{
public:
  /**
   * @brief		The constructor of strip
   * @param[in]	row : The rows of the strip
   * @param[in]	col : The columns of the strip
   * @param[in]	periodic : judge whether the boundray is periodic or not
   * @param[in]	Dirichlet : judge whether the boundray is Dirichlet or not
   * @param[in]	Gradient : judge whether the boundray is Neumann or not
   * @details
   * This is constructor. It is used to create an  object of the strip. Because
   * the strip have a boundary. So add 2 to the number of rows and columns passed in
   * @note	Notice that the value of row and col is bigger than 0.
   * @par Sample
   * @code
   *	strip A(10,10,false,true,false);
   *	cout<<A.row<<A.col<<endl;    // Both rows and columns of A are 12 
   * @endcode
   */
	strip(int row,int col, bool  periodic, bool  Dirichlet, bool Gradient):Matrix(row,col, periodic, Dirichlet,Gradient){}

	/**
    * @brief		The strip solver to deal with the wave equation using stripe decomposition
    * @param[in]	id: the id of the current process
    * @param[in]	p: The number of the whole processes
    * @param[in]	t_out: The time interval to generate the .dat file
    * @param[in]	t : The current time
    * @return       void
    * @attention	        
	* Notice that the number of the p is must smaller than the number of the strip's
	* rows, ie i_max. Because if p is bigger than the i_max, some process like process 0
	* has zero number of rows to deal with, which will lead to the wrong results.
    */
	void strip_solver(int id, int p, double t_out,double t);

	/**
    * @brief		This function is used to do strip decomposition
    * @param[in]	p : The number of the whole processes
	* @return       vector<vector<int>>
    * @details
    * This function returns the vector that includes two vectors. The first vector records 
	* the row start number in  the whole matrix for every process and the second vector 
	* records the number of rows that every process needs to deal with. Because for the matrix,
	* it has top boundary and bottom boundary, in the function, it will minus 2 from the number
	* of the strip rows.
    */
	vector<vector<int>> setup_partition(int p);


	/**
	* @brief		The evolution of the strip matrix
	* @param[in]    new_str: the pointer of the new_str whose type is strip
	* @param[in]    str: the pointer of the str whose type is strip
	* @param[in]    old_str: the pointer of the old_str whose type is strip
	* @param[in]    id: the id of the current process
	* @param[in]    p: The number of the whole processes
	* @param[in]    t: The time
	* @return       void
	* @details
	* The evolution is according to the wave equation. After evolution, according to the id and p 
	* value, decided to assign values to which boundary. Then, according to the different boundary
	* conditions to assign different values to the boundary.
	*/
	void strip_evolution(strip*& new_str, strip*& str, strip*& old_str, int id, int p, double& t);


	/**
    * @brief		Write the value of strip to the .dat file
    * @param[in]    str: the strip object
    * @param[in]    id: the id of the current process
    * @param[in]    chunk_num: Actually, it is the number of this str's rows minus 2
    * @param[in]    cnt: The number of the generated file
	* @param[in]    p: The number of the whole processes
    * @return       void
    * @details
    * According to different id values, determine the start line and end line of strip that write to 
	* the .dat file. 
    */
	void strip_to_file(strip str, int id, int chunk_num, int cnt, int p);

	/**
    * @brief		built the dataType of the strip. These dataType include strip_top and strip_bottom
    * @param[in]    str_value: one dimension array. It will receive the value of the strip
    * @param[in]    m: The number of the strip object's rows
    * @param[in]    n: The number of the strip object's cols
    * @return       void
    * @details
    * This function is used to create the datatype that used in MPI_Isend and MPI_Irecv. The strip_top only
	* sends row 1 of the strip. Because row 0 is used to receive the other stripe's boundary. The
	* strip_bottom only sends the second-to-last row of the strip, because the last row is used to receive the 
	* other strip's boundary
    */
	void built_stripDatatypes(double* str_value, int m, int n);
};


/**
* @class <Grid>
* @brief
* The class Grid is the subclass of the Matrix. The object of this class can solve wave function using
* grid decomposition.
*/
class Grid :public Matrix 
{

public:
	/**
    * @brief		The constructor of Grid
    * @param[in]	row : The rows of the Grid
    * @param[in]	col : The columns of the Grid
    * @param[in]	periodic : judge whether the boundray is periodic or not
    * @param[in]	Dirichlet : judge whether the boundray is Dirichlet or not
    * @param[in]	Gradient : judge whether the boundray is Neumann or not
    * @details
    * This is constructor. It is used to create an  object of the Grid. Because
    * the Grid have a boundary. So add 2 to the number of rows and columns passed in
    * @note	  Notice that the value of row and col is bigger than 0.
    * @par Sample
    * @code
    *	Grid A(10,10,false,true,false);
    *	cout<<A.row<<A.col<<endl;   // Both rows and columns of A are 12 
    * @endcode
    */
	Grid(int row, int col, bool  periodic, bool  Dirichlet, bool Gradient) :Matrix(row, col, periodic, Dirichlet, Gradient) {}


    /**
    * @brief		The grid solver to deal with the wave equation using grid decomposition
    * @param[in]	id: the id of the current process
    * @param[in]	p: the number of the whole processes
    * @param[in]	t_out: the time interval to generate the .dat file
    * @param[in]	t : The current time
    * @return       void
    * @attention
    * Notice that if the p is big enough, it will lead to some process deal with zero number of
    * rows or zero number cols, which will lead to the wrong results. So don't set p too large
    */
	void grid_solver(int id, int p, double t_out, double t);


    /**
    * @brief		This function is used to do grid decomposition
    * @param[in]	p_rows : the total number of processes in row direction
    * @param[in]	p_cols: the total number of processes in column direction
    * @return       vector<vector<int>>
    * @details
    * This function returns the vector that includes four vectors. The first vector records
    * the row start number in the whole matrix for every process and the second vector
    * records the number of rows that every process needs to deal with. The third vector records
    * the column start number in the whole matrix for every process and the fourth vector records
    * the number of columns that every process needs to deal with.  Because for the matrix,
    * it has top boundary and bottom boundary, in the function, it will minus 2 from the number of
    * the grid rows. Similarly, it will minus 2 from the number of the grid columns.
    */
	vector<vector<int>> setup_partition_grid(int p_rows, int p_cols);

	/**
    * @brief		The evolution of the Grid matrix
    * @param[in]    new_grid: the pointer of the new_grid whose type is Grid
    * @param[in]    grid: the pointer of the grid whose type is Grid
    * @param[in]    old_grid: the pointer of the old_grid whose type is Grid
    * @param[in]    row_index: the index of this id in row direction
    * @param[in]    col_index: the index of this id in column direction
    * @param[in]    p_rows: the total number of processes in row direction
    * @param[in]    p_cols: the total number of processes in column direction
    * @param[in]    t: the time
    * @return       void
    * @details
    * The evolution is according to the wave equation. After evolution, according to the row_index and col_index
    * value, decided to assign values to which boundary. Then, according to the different boundary conditions
    * to assign different values to the boundary.
    */
	void grid_evolution(Grid*& new_grid, Grid*& grid, Grid*& old_grid, int row_index, int col_index, int p_rows, int p_cols, double& t);

    /**
    * @brief		built the dataType of the grid. These dataType include grid_top, grid_bottom, grid_left, grid_right
    * @param[in]    grid_value: one dimension array. It will receive the value of the Grid
    * @param[in]    m: The number of the Grid object's rows
    * @param[in]    n: The number of the Grid object's cols
    * @return       void
    * @details
    * This function is used to create the datatype that used in MPI_Isend and MPI_Irecv. The grid_top only
    * sends row 1 of the Grid. Because row 0 is used to receive the other Grid's bottom boundary. The
    * grid_bottom only sends the second-to-last row of the Grid, because the last row is used to receive the
    * other Grid's top boundary. The grid_left only sends column 1 of the Grid, Because column 0 is used
    * to receive the other Grid's right boundary. The grid_right only sends the second-to-last column of the Grid, 
    * because the last column 0 is used to receive other Grid's left boundary
    */
	void built_gridDatatypes(double* grid_value, int m, int n);


    /**
    * @brief		Divide the total process p into p_rows×p_cols
    * @param[in]    p: The number of the whole processes
    * @param[in]    id: the id of the current process
    * @param[out]    p_rows: the total number of processes in row direction
    * @param[out]    p_cols: the total number of processes in column direction
    * @return       void
    */
	void find_dimensions(int p, int id, int& p_rows, int& p_cols);


    /**
    * @brief		Write the value of Grid to the .dat file
    * @param[in]    grid: the Grid object
    * @param[in]    p_rows: the total number of processes in row direction
    * @param[in]    p_cols: the total number of processes in column direction
    * @param[in]    row_index: the index of this id in row direction
    * @param[in]    col_index: the index of this id in column direction
    * @param[in]    cnt: The number of the generated file
    * @param[in]    id: the id of the current process
    * @return       void
    * @details
    * According to different id values, determine the start line and end line of Grid that write to
    * the .dat file.
    */
	void grid_to_file(Grid grid, int p_rows, int p_cols, int row_index, int col_index, int cnt, int id);
};
