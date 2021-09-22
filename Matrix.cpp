#include "Matrix.h"
#include <mpi.h>



// This function is used to initialize the matrix value
void Matrix::initial(int rows,int cols) {
	// set all values in the matrix to 0
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			this->values[i * cols + j] = 0;


	// sets half sinusoidal intitial disturbance - this is brute force - it can be done more elegantly
	for (int i = 1; i < rows - 1; i++) {
		for (int j = 1; j < cols - 1; j++) {
			double x = dx * i;
			double y = dy * j;
			double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

			if (dist < r_splash) {
				double h = 5.0 * (cos(dist / r_splash * PI) + 1.0);
				this->values[i * cols + j] = h;
			}
		}
	}
}

// built the Dirichlet boundray of the matrix
void Matrix::built_Dirichlet(int rows,int cols) {
	// Set the value of the 0th row and the last row of the matrix to 0
	for (int j = 0; j < cols; j++) {
		this->values[0 * cols + j] = 0;
		this->values[(rows - 1) * cols + j] = 0;
	}

	// Set the value of the 0th column and the last column of the matrix to 0
	for (int i = 1; i < rows - 1; i++) {
		this->values[i * cols] = 0;
		this->values[i * cols + cols - 1] = 0;
	}
}

// built the Neumann boundray of the matrix
void Matrix::built_Gradient(int rows, int cols) {
	// The value of column 0 of the matrix is equal to the value of column 1 of the matrix
	// The value of the last column of the matrix is equal to the value of the second-to-last
	// column of the matrix
	for (int i = 0; i < rows; i++) {
		this->values[i * cols] = this->values[i * cols + 1];
		this->values[i * cols + cols - 1] = this->values [i * cols + cols - 2];
	}

	// The value of row 0 of the matrix is equal to the value of row 1 of the matrix
	// The value of the last row of the matrix is equal to the value of the second-to-last
	// row of the matrix
	for (int j = 0; j < cols; j++)
	{
		this->values[j] = this->values[1 * cols + j];
		this->values[(rows - 1) * cols + j] = this->values[(rows - 2) * cols + j];
	}
}

// built the periodic boundray of the matrix
void Matrix::built_Periodic(int rows, int cols) {
	// The value of the row 0 of the matrix is equal to the value of the second-to-last
	// row of the matrix. The value of the last row of the matrix is equal to the value 
	// of the row 1 of the matrix
	for (int j = 1; j < cols - 1; j++) {
		this->values[0 * cols + j] = this->values[(rows - 2) * cols + j]; 
		this->values[(rows - 1) * cols + j] = this->values[1 * cols + j];
	}

	// The value of column 0 of the matrix is equal to the value of second-to-last column of the matrix
	// The value of the last column of the matrix is equal to the value of the column 1 of the matrix
	for (int i = 1; i < rows - 1; i++) {
		this->values[i * cols + 0] = this->values[i * cols + cols - 2];
		this->values[i * cols + cols - 1] = this->values[i * cols + 1];
	}

	// The calculation does not require the four vertices of the matrix. 
	// So set the four vertices of the matrix to 0
	this->values[0] = 0;
	this->values[0 * cols + cols - 1] = 0;
	this->values[(rows - 1) * cols + 0] = 0;
	this->values[(rows - 1) * cols + cols - 1] = 0;

}



Matrix::Matrix(int row, int col, bool  periodic, bool  Dirichlet,bool Gradient) {
	this->row = row + 2;
	this->col = col + 2;
	// Set the boundary conditions. We only allow one boundary condition
	// If all three boundary conditions are false, we default to Dirichlet boundary
	if (periodic) {
		this->periodic = true;
		this->Dirichlet = false;
		this->Gradient = false;
	}else if (Gradient) {
		this->Gradient = true;
		this->periodic = false;
		this->Dirichlet = false;
	}
	else {
		this->Dirichlet = true;
		this->periodic = false;
		this->Gradient = false;
	}
	// Create a one-dimensional array. The size is rows multiply columns
	this->values = new double[this->row * this->col];
	// Initialize the value of the matrix
	initial(this->row, this->col);
}



// Print the valiue of the matrix
void Matrix::printMatrix()
{
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->col; j++)
		{
			cout << this->values[i * this->col + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// Acoording to the wave equation, do the evolution
void Matrix::evolution( Matrix*& new_grid,Matrix*& grid, Matrix*& old_grid,int id,double& t) {
	for (int i = 1; i < grid->row - 1; i++) {
		for (int j = 1; j < grid->col - 1; j++) {
			new_grid->values[i * grid->col + j] = (pow(dt * c, 2.0) * ((grid->values[(i + 1) * grid->col + j] -
				2.0 * grid->values[i * grid->col + j] + grid->values[(i - 1) * grid->col + j]) / pow(dx, 2.0) +
				(grid->values[i * grid->col + j + 1] - 2.0 * grid->values[i * grid->col + j] +
				grid->values[i * grid->col + j - 1]) / pow(dy, 2.0)) + 2.0 * grid->values[i * grid->col + j]
			    - old_grid->values[i * grid->col + j]);
		}
	}

	// Then assign boundary values according to boundary conditions
	if (this->periodic) {
		new_grid->built_Periodic(this->row, this->col);
	}else if (this->Gradient){
		new_grid->built_Gradient(this->row, this->col);
	}
	else {
		new_grid->built_Dirichlet(this->row, this->col);
	}

	// record the pointer of the old_grid
	Matrix* old_address = old_grid;

	// Assign the pointer value of the grid to old_grid
	old_grid = grid;
	// Assign the pointer value of the new_grid to grid
	grid = new_grid;
	// Assign the address of the original old_grid to new_grid
	new_grid = old_address;

	// After change the pointer, the grid have the new_grid value , the old_grid has the grid value. 
	// the new_grid value will be updated next iteration

}