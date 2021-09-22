#include "Solver.h"
#include <mpi.h>


MPI_Datatype strip_top,strip_bottom;
MPI_Datatype grid_top, grid_bottom, grid_left, grid_right;

// built the dataType of the strip. These dataType include strip_top and strip_bottom
void strip::built_stripDatatypes(double* str_value, int m, int n) {
	MPI_Aint add_start;

	int block_length = n; // The number of elements to send
	MPI_Datatype typevalue = MPI_DOUBLE;  // The type
	MPI_Aint address;

	MPI_Get_address(str_value, &add_start);
	// Because row 0 is added later, the top data to send in row 1
	MPI_Get_address(&str_value[1 * n], &address);
	address = address - add_start;
	// create the strip_top dataType and commit
	MPI_Type_create_struct(1, &block_length, &address, &typevalue, &strip_top);
	MPI_Type_commit(&strip_top);

	// Because last row is added later, the bottom data to send in the second-to-last row 
	MPI_Get_address(&str_value[(m - 2) * n], &address);
	address = address - add_start;
	// create the strip_bottom dataType and commit
	MPI_Type_create_struct(1, &block_length, &address, &typevalue, &strip_bottom);
	MPI_Type_commit(&strip_bottom);
}

// Write the value of strip to the .dat file
void strip:: strip_to_file(strip str,int id, int chunk_num, int cnt,int p) {
	stringstream fname;
	fstream f1;
	// The path
	fname << "../../out/output" << "_" << cnt <<"th_"<<id<<"of"<<p<<"_"<<imax<<"x"<<jmax<< ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	if (id == 0) {
		// Start print from row 0 to second-to-last row
		for (int i = 0; i < chunk_num + 1; i++) {
			for (int j = 0; j < str.col; j++) {
				f1 << str.values[i * str.col + j] << "\t";
			}
			f1 << endl;
		}
		f1.close();
	}else if (id == p - 1) {
		// Start print from row 1 to last row
		for (int i = 1; i < chunk_num + 2; i++) {
			for (int j = 0; j < str.col; j++) {
				f1 << str.values[i * str.col + j] << "\t";
			}
			f1 << endl;
		}
		f1.close();
	}
	else {
		// Start print from row 1 to second-to-last row
		for (int i = 1; i < chunk_num + 1; i++) {
			for (int j = 0; j < str.col; j++) {
				f1 << str.values[i * str.col + j] << "\t";
			}
			f1 << endl;
		}
		f1.close();
	}
}

// This function is used to do strip decomposition for the matrix
vector<vector<int>> strip::setup_partition(int p) {
	vector<vector<int>> result;
	vector<int> pro_chunk;
	vector<int> row_star;
	// This array is used to store row start number in whole matrix for every process
	int* row_start_list = new int[p];
	// This array is used to store the number of rows that every process need to deal with
	int* process_chunk_list = new int[p];


	// Becasue the first row and last row are boundray, the number of rows to deal with
	// for all processes are the total rows minus 2
	int rows_left = this->row - 2;
	row_start_list[0] = 0;
	for (int n = 0; n < p - 1; n++)
	{
		int rows_assigned = rows_left / (p - n);
		rows_left -= rows_assigned;
		row_start_list[n + 1] = row_start_list[n] + rows_assigned;
		process_chunk_list[n] = rows_assigned;
	}
	// Hand over all the remaining rows to the last process
	process_chunk_list[p - 1] = this->row - 2 - row_start_list[p - 1];

	// Pass the data of the array to the vector
	for (int i = 0; i < p; i++) {
		pro_chunk.push_back(process_chunk_list[i]);
	}

	for (int i = 0; i < p; i++) {
		row_star.push_back(row_start_list[i]);
	}
	
	result.push_back(row_star);
	result.push_back(pro_chunk);

	// release the memory
	delete[] row_start_list;
	delete[] process_chunk_list;
	row_start_list = nullptr;
	process_chunk_list=nullptr;
	return result;
}

// The strip solver to deal with the wave equation using strip decomposition
void strip::strip_solver(int id, int p, double t_out, double t) {

	vector<vector<int>> result = setup_partition(p);
	// This vector stores row start number in whole matrix for every process
	vector<int> row_start = result[0];
	// This vector stores the number of rows that every process need to deal with
	vector<int> process_chunk = result[1];
	
	// Print the warning and decomposition information
	if (id == 0) {
		if (p > imax - 2) {
			cout << "The input processes are bigger than the rows, the strip decomposition will lead to wrong results" << endl;
			cout << "Please input the process that smaller than the rows to get the correct answer" << endl;
		}
		for (int i = 0; i < p; i++) {
			cout << "p" << i << ": " << row_start[i] << " " << process_chunk[i] << endl;
		}
	}

	// Create a matrix whose rows and cols are same to this matrix. Because the rows and cols are input
	// numbers added 2, here should be minus 2. This object is used to initialize the value of strip
	// type object
	Matrix A(this->row - 2, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);
	
	// Because the first row is boundray, the actual row start num in the whole matrix of 
	// this process should be added 1
	int row_stnum = row_start[id] + 1;
	// This parameter store the number of rows that this process will deal with
	int chunk_num = process_chunk[id];
	// create three strip objects. Using pointer is easy to swap values later
	strip* new_stri = new strip(chunk_num, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);
	strip* stri = new strip(chunk_num, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);
	strip* old_stri = new strip(chunk_num, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);
	
	// Initial the value of the stri and old_stri. Because the first and last rows 
	// of stri and old_stri are used to receive boundaries, assign values from the 
	// second row to second-to-last row. 
	for (int i = row_stnum; i < row_stnum + chunk_num; i++) {
		for (int j = 0; j < this->col; j++) {
			stri->values[(i - row_stnum + 1) * stri->col + j] = A.values[i * A.col + j];
			old_stri->values[(i - row_stnum + 1) * stri->col + j] = A.values[i * A.col + j];
		}
	}

	// the id number of the last id
	int idOfup;
	// the id number of the next id
	int idOfdown;

	// general condition
	idOfup = id - 1;
	idOfdown = id + 1;

	
	if (id == 0) {
		// if periodic, the idofup of process 0 is equal to the last process
		if (this->periodic) {
			idOfup = p - 1;
		}
		else {
			// if not periodic, the process will not receive any message, set idofup is MPI_PROC_NULL
			// this means empty process. Use empty process communication to do nothing
			idOfup = MPI_PROC_NULL;
		}
	}

	if (id == p - 1) {
		// if periodic, the idofdown of last process is equal to process 0
		if (this->periodic) {
			idOfdown = 0;
		}
		else {
			idOfdown = MPI_PROC_NULL;
		}
	}

	// genarte the first .dat file to record the initial state before the evolution
	int out_cnt = 0, it = 0;
	this->strip_to_file(*stri, id, chunk_num, out_cnt, p);
	out_cnt++;
	t_out += dt_out;

	MPI_Request* request = new MPI_Request[8];
	
	while (t < t_max) 
	{
		// create our own datatype
		built_stripDatatypes(stri->values, stri->row, stri->col);
		built_stripDatatypes(old_stri->values, old_stri->row, old_stri->col);

		// send the boundary value to corresponding process
		MPI_Isend(stri->values, 1, strip_top, idOfup, 0, MPI_COMM_WORLD, &request[0]);
		MPI_Isend(stri->values, 1, strip_bottom, idOfdown, 1, MPI_COMM_WORLD, &request[1]);
		MPI_Isend(old_stri->values, 1, strip_top, idOfup, 2, MPI_COMM_WORLD, &request[4]);
		MPI_Isend(old_stri->values, 1, strip_bottom, idOfdown, 3, MPI_COMM_WORLD, &request[5]);

		// Create a one-dimensional array for receiving data
		double* recv_top_stri;  // the strip_bottom type. The second-to-last row has values
		double* recv_bottom_stri;  // the strip_top type. The row 1 has values
		double* recv_top_oldstri;
		double* recv_bottom_oldstri;
		if (id == 0) {
			// if periodic, it will receive the value from the last process and the size 
			// of the array is temp_rows * stri->col. if not periodic, it won't receive 
			// the values. we also set this size, which has no effect on the result
			int temp_rows = max<int>(process_chunk[p - 1], process_chunk[0]) + 2;
			recv_top_stri = new double[temp_rows * stri->col];
			recv_top_oldstri = new double[temp_rows * stri->col];
		}
		else {
			// if the process is not zero, we choose the max function to avoid overflow
			int temp_rows = max<int>(process_chunk[id - 1], process_chunk[id]) + 2;
			recv_top_stri = new double[temp_rows * stri->col];
			recv_top_oldstri = new double[temp_rows * stri->col];
		}

		if (id == p - 1) {
			// if periodic, it will receive the value from the process 0 and the size 
			// of the array is temp_rows * stri->col. if not periodic, it won't receive 
			// the values. we also set this size, which has no effect on the result
			int temp_rows = max<int>(process_chunk[0], process_chunk[p - 1]) + 2;
			recv_bottom_stri = new double[temp_rows * stri->col];
			recv_bottom_oldstri = new double[temp_rows * stri->col];
		}
		else {
			// if the process is not the last process, we choose the max function to avoid overflow
			int temp_rows = max<int>(process_chunk[id + 1], process_chunk[id]) + 2;
			recv_bottom_stri = new double[temp_rows * stri->col];
			recv_bottom_oldstri = new double[temp_rows * stri->col];
		}
		
	
		MPI_Irecv(recv_top_stri, 1, strip_bottom, idOfup, 1, MPI_COMM_WORLD, &request[2]);
		MPI_Irecv(recv_bottom_stri, 1, strip_top, idOfdown, 0, MPI_COMM_WORLD, &request[3]);
		MPI_Irecv(recv_top_oldstri, 1, strip_bottom, idOfup, 3, MPI_COMM_WORLD, &request[6]);
		MPI_Irecv(recv_bottom_oldstri, 1, strip_top, idOfdown, 2, MPI_COMM_WORLD, &request[7]);

		MPI_Waitall(8, request, MPI_STATUS_IGNORE);
		
		 
		// Assign the row 1 values of recv_ bottom and recv_bottom_oldstri to the last row of 
		// the stri and old_stri
		if (this->periodic) {
			for (int j = 0; j < stri->col; j++) {
				stri->values[(stri->row - 1) * stri->col + j] = recv_bottom_stri[1 * stri->col + j];
				old_stri->values[(old_stri->row - 1) * old_stri->col + j] = recv_bottom_oldstri[1 * old_stri->col + j];
			}
		}
		else {
			// if not periodic and process p - 1, the recv_bottom_stri and oldstri won't receive the values
			// from other process.
			if (id != p - 1) {
				for (int j = 0; j < stri->col; j++) {
					stri->values[(stri->row - 1) * stri->col + j] = recv_bottom_stri[1 * stri->col + j];
					old_stri->values[(old_stri->row - 1) * old_stri->col + j] = recv_bottom_oldstri[1 * old_stri->col + j];
				}
			}
		}
	
		// Assign the value of the second-to-last row of recv_top_stri and oldstri to the row 0 of stri
		// and old_stri
		if (this->periodic) {
			int temp_rows = process_chunk[id] + 2;
			if (id == 0) {
				//temp_rows = max<int>(process_chunk[p - 1], process_chunk[0]) + 2;
				for (int j = 0; j < stri->col; j++) {
					stri->values[j] = recv_top_stri[(temp_rows - 2) * stri->col + j];
					old_stri->values[j] = recv_top_oldstri[(temp_rows - 2) * old_stri->col + j];
				}
			}
			else {
				//temp_rows = max<int>(process_chunk[id - 1], process_chunk[id]) + 2;
				for (int j = 0; j < stri->col; j++) {
					stri->values[j] = recv_top_stri[(temp_rows - 2) * stri->col + j];
					old_stri->values[j] = recv_top_oldstri[(temp_rows - 2) * old_stri->col + j];
				}
			}
		}
		else {
			// For process 0, if not periodic, the recv_top_stri and old_stri won't receive 
			// the value from other processes
			if (id != 0) {
				int temp_rows = process_chunk[id] + 2;
				for (int j = 0; j < stri->col; j++) {
					stri->values[j] = recv_top_stri[(temp_rows - 2) * stri->col + j];
					old_stri->values[j] = recv_top_oldstri[(temp_rows - 2) * old_stri->col + j];
				}
			}
		}
		// do evolution
		this->strip_evolution(new_stri, stri, old_stri, id, p,t);
	    // write the value to .dat file
		if (t_out <= t) {
			 this->strip_to_file(*stri, id, chunk_num, out_cnt, p);
			out_cnt++;
			t_out += dt_out;
		}
		it++;
		
		// release memory
		delete[] recv_top_stri;
		delete[] recv_bottom_stri;
		delete[] recv_top_oldstri;
		delete[] recv_bottom_oldstri;
		recv_top_stri = nullptr;
		recv_bottom_stri = nullptr;
		recv_top_oldstri = nullptr;
		recv_bottom_oldstri = nullptr;

		MPI_Type_free(&strip_top);
		MPI_Type_free(&strip_bottom);

		// add Barrier to wait all processes. Then do the next evolution
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (id == 0) {
		cout << "output: " << out_cnt - 1 << "\tt: " << t << "\titeration: " << it - 1 << endl;
	}

	// release memory
	delete[] request;
	delete new_stri;
	delete stri;
	delete old_stri;
	request = nullptr;
	new_stri = nullptr;
	stri = nullptr;
	old_stri = nullptr;
}



// The evolution is according to the wave equation
void strip::strip_evolution(strip*& new_str, strip*& str, strip*& old_str,int id,int p,double& t) {
	for (int i = 1; i < new_str->row - 1; i++)
		for (int j = 1; j < new_str->col - 1; j++) {
			new_str->values[i*str->col+j] = pow(dt*c,2.0) * ((str->values[(i + 1)*str->col+j] -
				2.0 * str->values[i*str->col+j] + str->values[(i - 1)*str->col+j]) / pow(dx, 2.0) +
				(str->values[i*str->col+j + 1] - 2.0 * str->values[i*str->col+j] +
				str->values[i*str->col+j - 1]) / pow(dy, 2.0)) + 2.0 * str->values[i*str->col+j]
				- old_str->values[i*str->col+j];
			
		}

	// Determine the distribution method of boundary values according to boundary conditions
	if (this->periodic) {
		// Assign the left and right columns values
		for (int i = 0;i < str->row; i++) {
			// left column
			new_str->values[i * str->col] = new_str->values[i * str->col + str->col - 2];
			// right column
			new_str->values[i * str->col + str->col - 1] = new_str->values[i * str->col + 1];

		}
		// if id=0 or p-1, the top and bottom value will receive from other process. Therefore, we
		// won't assign values to the top and bottom
	}
	else if (this->Gradient) {
		// Assign the left and right columns values
		for (int i = 0; i < str->row; i++) {
			// left column
			new_str->values[i * str->col + 0] = new_str->values[i * str->col + 1];
            // right column
			new_str->values[i * str->col + str->col - 1] = new_str->values[i * str->col + str->col - 2];
		}
		// if id = 0, we assign the top value of the matrix
		if (id == 0) {
			for (int j = 0; j < str->col; j++) {
				new_str->values[j] = new_str->values[1 * str->col + j];
			}
		}
		// if id = p - 1, we assign the bottom value of the matrix
		if (id == p - 1) {
			for (int j = 0; j < str->col; j++) {
				new_str->values[(str->row - 1) * str->col + j] = new_str->values[(str->row - 2) * str->col + j];
			}
		}
	} else {
		// The default boundary condition is  Dirichlet
		// Assign values of the left and right colunmns
		for (int i = 0; i< str->row; i++) {
			// left colunmn
			new_str->values[i * str->col] = 0;
            // right colunmn
			new_str->values[i * str->col + str->col - 1] = 0;
		
		}
		// if id = 0, we assign the top value of the matrix
		if (id == 0) {
			for (int j = 0; j < str->col; j++) {
				new_str->values[j] = 0;
			}
		} // if id = p - 1, we assign the bottom value of the matrix
		else if (id == p - 1) {
			for (int j = 0; j < str->col; j++) {
				new_str->values[(str->row - 1) * str->col + j] = 0;
			}
		}
	}
	t += dt;
	// swap the address to make the stri has the updated value
	strip* old_address = old_str;
	old_str = str;
	str = new_str;
	new_str = old_address;
}

// built the dataType of the grid. These dataType include grid_top, grid_bottom, grid_left, grid_right
void Grid::built_gridDatatypes(double* grid_value, int m, int n) {
	vector<int> block_lengths;
	vector<MPI_Datatype> typelist;
	vector<MPI_Aint> addresses;
	MPI_Aint add_start;

	// left type
	// Because column 0 is added later, the left data to send in column 1
	for (int i = 0; i < m; i++) {
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&grid_value[i * n + 1], &temp_address);
		addresses.push_back(temp_address);
	}
	MPI_Get_address(grid_value, &add_start);
	for (int i = 0; i < m; i++)
		addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &grid_left);
	MPI_Type_commit(&grid_left);

	// right type
	// Because right column is added later, the right data to send in second-to-last column
	block_lengths.resize(0);
	typelist.resize(0);
	addresses.resize(0);
	for (int i = 0; i < m; i++) {
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&grid_value[i * n + n - 2], &temp_address);
		addresses.push_back(temp_address);
	}
	for (int i = 0; i < m; i++)
		addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &grid_right);
	MPI_Type_commit(&grid_right);

	// top type
	int block_length = n;
	MPI_Datatype typevalue = MPI_DOUBLE;
	MPI_Aint address;
	// Because row 0 is added later, the top data to send in row 1
	MPI_Get_address(&grid_value[1 * n], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typevalue, &grid_top);
	MPI_Type_commit(&grid_top);

	// bottom type
	// Because last row is added later, the bottom data to send in the second-to-last row 
	MPI_Get_address(&grid_value[(m - 2) * n], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typevalue, &grid_bottom);
	MPI_Type_commit(&grid_bottom);
	
}

// Write the value of Grid to the .dat file
void Grid::grid_to_file(Grid grid, int p_rows, int p_cols,int row_index, int col_index, int cnt,int id){
	stringstream fname;
	fstream f1;
	// path
	fname << "../../out/output" << "_" << cnt << "th_" << id << "of(" <<row_index<<","<<col_index << ")_" << imax << "x" << jmax << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	if (row_index == 0 && col_index == 0) {
		// upper left process
		// rows from 0 to the second-to-last row
		// columns from 0 to the second-to-lat column
		for (int i = 0; i < grid.row - 1; i++) {
			for (int j = 0; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (row_index == 0 && col_index == p_cols - 1) {
		// upper right process
		// rows from 0 to the second-to-last row
		// columns from 1 to the last column
		for (int i = 0; i < grid.row - 1; i++) {
			for (int j = 1; j < grid.col; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (row_index == p_rows - 1 && col_index == 0) {
		// bottom left process
		// rows from 1 to the last
		// columns from 0 to the second-to-last
		for (int i = 1; i < grid.row; i++) {
			for (int j = 0; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (row_index == p_rows - 1 && col_index == p_cols - 1) {
		// bottom right process
		// rows from 1 to the second-to-last row
		// columns from 1 to the last column
		for (int i = 1; i < grid.row; i++) {
			for (int j = 1; j < grid.col; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (row_index == 0) {
		// the top processes
		// rows from 0 to the second-to-last
		// columns from 1 to second-to-last
		for (int i = 0; i < grid.row - 1; i++) {
			for (int j = 1; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (row_index == p_rows - 1) {
		// the bottom processes
		// rows from 1 to the last
		// columns from 1 to the second-to-last
		for (int i = 1; i < grid.row; i++) {
			for (int j = 1; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (col_index == 0) {
		// the left processes
		// rows from 1 to second-to-last
		// columns from 0 to second-to-last
		for (int i = 1; i < grid.row - 1; i++) {
			for (int j = 0; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else if (col_index == p_cols - 1) {
		// the right processes
		// rows from 1 to the second-to-last
		// columns from 1 to the last
		for (int i = 1; i < grid.row - 1; i++) {
			for (int j = 1; j < grid.col; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	else {
		// the general processes
		// rows from 1 to the second-to-last
		// columns from 1 to the second-to-last
		for (int i = 1; i < grid.row - 1; i++) {
			for (int j = 1; j < grid.col - 1; j++) {
				f1 << grid.values[i * grid.col + j] << "\t";
			}
			f1 << endl;
		}
	}
	
	f1.close();
}

// Divide the total process p into p_rows×p_cols
void Grid::find_dimensions(int p, int id, int& p_rows, int& p_cols) {
	int min_gap = p;
	int top = sqrt(p) + 1;
	for (int i = 1; i <= top; i++)
	{
		if (p % i == 0)
		{
			int gap = abs(p / i - i);
			if (gap < min_gap)
			{
				min_gap = gap;
				p_rows = i;
				p_cols = p / i;
			}
		}
	}

	if (id == 0)
		cout << "Divide " << p << " into " << p_rows << " by " << p_cols << " grid" << endl;
}

// This function is used to do grid decomposition
vector<vector<int>> Grid::setup_partition_grid(int p_rows,int p_cols) {
	vector<vector<int>> result;
	vector<int> row_start;
	vector<int> col_start;
	vector<int> row_process_chunk;
	vector<int> col_process_chunk;

	// This array is used to store row start number in whole matrix for every process
	int* row_start_list = new int[p_rows];
	// This array is used to store the number of rows that every process need to deal with
	int* row_process_chunk_list = new int[p_rows];

	// Becasue the first row and last row are boundray, the number of rows to deal with
    // for all processes are the total rows minus 2
	int rows_left = this->row - 2;
	row_start_list[0] = 0;
	for (int n = 0; n < p_rows - 1; n++)
	{
		int rows_assigned = rows_left / (p_rows - n);
		rows_left -= rows_assigned;
		row_start_list[n + 1] = row_start_list[n] + rows_assigned;
		row_process_chunk_list[n] = rows_assigned;
	}
    // Hand over all the remaining rows to the last p_col processes
	row_process_chunk_list[p_rows - 1] = this->row - 2 - row_start_list[p_rows - 1];

	// Pass the data of the array to the vector
	for (int i = 0; i < p_rows; i++) {
		row_process_chunk.push_back(row_process_chunk_list[i]);
	}

	for (int i = 0; i < p_rows; i++) {
		row_start.push_back(row_start_list[i]);
	}

	result.push_back(row_start);
	result.push_back(row_process_chunk);

	// This array is used to store column start number in whole matrix for every process
	int* col_start_list = new int[p_cols];
	// This array is used to store the number of columns that every process need to deal with
	int* col_process_chunk_list = new int[p_cols];

	// Becasue the first column and last column are boundray, the number of columns to deal with
	// for all processes are the total columns minus 2
	int cols_left = this->col - 2;
	col_start_list[0] = 0;
	for (int n = 0; n < p_cols - 1; n++)
	{
		int cols_assigned = cols_left / (p_cols - n);
		cols_left -= cols_assigned;
		col_start_list[n + 1] = col_start_list[n] + cols_assigned;
		col_process_chunk_list[n] = cols_assigned;
	}
	// Hand over all the remaining columns to the last p_rows processes
	col_process_chunk_list[p_cols - 1] = this->col - 2 - col_start_list[p_cols - 1];

	// Pass the data of the array to the vector
	for (int i = 0; i < p_cols; i++) {
		col_process_chunk.push_back(col_process_chunk_list[i]);
	}

	for (int i = 0; i < p_cols; i++) {
		col_start.push_back(col_start_list[i]);
	}
	result.push_back(col_start);
	result.push_back(col_process_chunk);

	// release the memory
	delete[] row_start_list;
	delete[] row_process_chunk_list;
	delete[] col_start_list;
	delete[] col_process_chunk_list;
	row_start_list = nullptr;
	row_process_chunk_list = nullptr;
	col_start_list = nullptr;
	col_process_chunk_list = nullptr;
	return result;

}

// The  grid solver to deal with the wave equation using grid decomposition
void Grid::grid_solver(int id, int p, double t_out, double t) {
	// the total number of processes in row direction
	int p_rows;
	// the total number of processes in column direction
	int p_cols;
	// This function is used to get the p_rows and p_cols value
	this->find_dimensions(p, id, p_rows, p_cols);

	// If p_rows or p_cols equalling 1, call strip decomposition
	if (p_rows == 1 || p_cols == 1) {
		strip S(this->row - 2, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);
		S.strip_solver(id, p, t, t_out);
		return;
	}

	vector<vector<int>> result = this->setup_partition_grid(p_rows,p_cols);
	// This vector stores row start number in whole matrix for every process
	vector<int> row_start = result[0];
	// This vector stores the number of rows that every process need to deal with
	vector<int> row_process_chunk = result[1];
	// This vector stores column start number in whole matrix for every process
	vector<int> col_start = result[2];
	// This vector stores the number of columns that every process need to deal with
	vector<int> col_process_chunk = result[3];

	// The index of the id in row dirction
	int row_index = id / p_cols;
	// The index of the id in column dirction
	int col_index = id % p_cols;

	// Print the warning and decomposition information
    if (id == 0) {
		if (row_process_chunk[0] == 0 || col_process_chunk[0] == 0) {
			cout << "The number of processes is too bigger, it will lead to the wrong result." << endl;
			cout << "Please input the smaller number of processes to get the correct answer" << endl;
		}
		for (int i = 0; i < p; i++) {
			cout << "p: " << i << "  row_start " << row_start[i / p_cols] <<  "  col start " << col_start[i % p_cols] << " row_chunk: " << row_process_chunk[i / p_cols] << " col_chunk: " << col_process_chunk[i % p_cols] << endl;
		}
	}


	// Create a matrix whose rows and cols are same to this matrix. Because the rows and cols are input
	// numbers added 2, here should be minus 2. This object is used to initialize the value of strip
	// type object
	Matrix A(this->row - 2, this->col - 2, this->periodic, this->Dirichlet, this->Gradient);

	// Because the first row is boundray, the actual row start num in the whole matrix of 
    // this process should be added 1. The same to the column start number.
	int row_stnum = row_start[row_index] + 1;
	int col_stnum = col_start[col_index] + 1;
	// These parameters store the number of rows and columns that this process will deal with
	int row_chunk_num = row_process_chunk[row_index];
	int col_chunk_num = col_process_chunk[col_index];
	// create three Grid objects. Using pointer is easy to swap values later
	Grid* new_grid = new Grid(row_chunk_num, col_chunk_num, this->periodic, this->Dirichlet, this->Gradient);
	Grid* grid = new Grid(row_chunk_num, col_chunk_num, this->periodic, this->Dirichlet, this->Gradient);
	Grid* old_grid = new Grid(row_chunk_num, col_chunk_num, this->periodic, this->Dirichlet, this->Gradient);

    // Initial the value of the grid and old_grid
	for (int i = row_stnum; i < row_stnum + row_chunk_num; i++) {
		for (int j = col_stnum; j < col_stnum + col_chunk_num; j++) {
			// 从第一行第一列开始
			grid->values[(i - row_stnum + 1) * grid->col + j - col_stnum + 1] = A.values[i * A.col + j];
			old_grid->values[(i - row_stnum + 1) * grid->col + j - col_stnum + 1] = A.values[i * A.col + j];
		}
	}
   

	// the id number of the up id
	int idOfup;
	// the id number of the down id
	int idOfdown;
	// the id number of the left id
	int idOfleft;
	// the id number of the right id
	int idOfright;

	// special condition for periodic boundary condition
	if (p == 1 && this->periodic)
	{
		idOfup = 0;
		idOfdown = 0;
		idOfleft = 0;
		idOfright = 0;
	}


	// general condition
	idOfup = id - p_cols;
	idOfdown = id + p_cols;
	idOfleft = id - 1;
	idOfright = id + 1;

	if (row_index == 0) {
		if (this->periodic) {
			idOfup = (p_rows - 1) * p_cols + id;
		}
		else {
			// if not periodic, the processes with row_index equalling 0 will not receive any message, 
			// set idofup is MPI_PROC_NULL, this means empty process. Use empty process communication to do nothing
			idOfup= MPI_PROC_NULL;
		}
	}

	if (row_index == p_rows - 1) {
		if (this->periodic) {
			idOfdown = id - (p_rows - 1) * p_cols;
		}
		else {
			idOfdown= MPI_PROC_NULL;
		}
	}

	if (col_index == 0) {
		if (this->periodic) {
			idOfleft = id + p_cols - 1;
		}
		else {
			idOfleft = MPI_PROC_NULL;
		}
	}

	if (col_index == p_cols - 1) {
		if (this->periodic) {
			idOfright = id - (p_cols - 1);
		}
		else {
			idOfright = MPI_PROC_NULL;
		}
	}
	
	// genarte the first .dat file to record the initial state before the evolution
	int out_cnt = 0, it = 0;
	this->grid_to_file(*grid,p_rows,p_cols, row_index, col_index, out_cnt, id);
	out_cnt++;
	t_out += dt_out;

	MPI_Request* request = new MPI_Request[16];
	while ( t < t_max) {
		// create our own datatype
		built_gridDatatypes(grid->values, grid->row, grid->col);
		built_gridDatatypes(old_grid->values, old_grid->row, old_grid->col);

		// send the boundary value to corresponding process
		MPI_Isend(grid->values, 1, grid_top, idOfup, 0, MPI_COMM_WORLD, &request[0]);
		MPI_Isend(grid->values, 1, grid_bottom, idOfdown, 1, MPI_COMM_WORLD, &request[1]);
		MPI_Isend(grid->values, 1, grid_left, idOfleft, 2, MPI_COMM_WORLD, &request[2]);
		MPI_Isend(grid->values, 1, grid_right, idOfright, 3, MPI_COMM_WORLD, &request[3]);

		MPI_Isend(old_grid->values, 1, grid_top, idOfup, 4, MPI_COMM_WORLD, &request[4]);
		MPI_Isend(old_grid->values, 1, grid_bottom, idOfdown, 5, MPI_COMM_WORLD, &request[5]);
		MPI_Isend(old_grid->values, 1, grid_left, idOfleft, 6, MPI_COMM_WORLD, &request[6]);
		MPI_Isend(old_grid->values, 1, grid_right, idOfright, 7, MPI_COMM_WORLD, &request[7]);

		// Create a one-dimensional array for receiving data
		double* recv_top_grid; 
		double* recv_bottom_grid;
		double* recv_left_grid;
		double* recv_right_grid;

		double* recv_top_oldgrid; 
		double* recv_bottom_oldgrid;
		double* recv_left_oldgrid;
		double* recv_right_oldgrid;

		if (row_index == 0) {
		    // if periodic, the processes with row_index equalling 0 will receive the value from the processes 
			// with row_index eqaulling  p_rows -1, and the size of the array is temp_rows * grid->col. if not
			// periodic, it won't receive the values. we also set this size, which has no effect on the result
			int temp_rows = max<int>(row_process_chunk[p_rows - 1], row_process_chunk[row_index]) + 2;
			recv_top_grid = new double[temp_rows * grid->col];
			recv_top_oldgrid = new double[temp_rows * grid->col];
		}
		else {
			// if the row_index of the process is not zero, we choose the max function to avoid overflow
			int temp_rows = max<int>(row_process_chunk[row_index - 1], row_process_chunk[row_index]) + 2;
			recv_top_grid = new double[temp_rows * grid->col];
			recv_top_oldgrid = new double[temp_rows * grid->col];
		}

		if (row_index == p_rows - 1) {
			// if periodic, the processes with row_index equalling p_rows - 1  will receive the value from the processes 
			// with row_index eqaulling  0, and the size of the array is temp_rows * grid->col. if not
			// periodic, it won't receive the values. we also set this size, which has no effect on the result
			int temp_rows = max<int>(row_process_chunk[0], row_process_chunk[row_index]) + 2;
			recv_bottom_grid = new double[temp_rows * grid->col];
			recv_bottom_oldgrid = new double[temp_rows * grid->col];
		}
		else {
			int temp_rows = max<int>(row_process_chunk[row_index + 1], row_process_chunk[row_index]) + 2;
			recv_bottom_grid = new double[temp_rows * grid->col];
			recv_bottom_oldgrid = new double[temp_rows * grid->col];
		}

		if (col_index == 0) {
		    // the same concept to rows
			int temp_cols = max<int>(col_process_chunk[p_cols - 1], col_process_chunk[col_index]) + 2;
			recv_left_grid = new double[grid->row * temp_cols];
			recv_left_oldgrid = new double[grid->row * temp_cols];
		}
		else {
			int temp_cols = max<int>(col_process_chunk[col_index - 1], col_process_chunk[col_index]) + 2;
			recv_left_grid = new double[grid->row * temp_cols];
			recv_left_oldgrid = new double[grid->row * temp_cols];
		}

		if (col_index == p_cols - 1) {
			int temp_cols = max<int>(col_process_chunk[0], col_process_chunk[col_index]) + 2;
			recv_right_grid = new double[grid->row * temp_cols];
			recv_right_oldgrid = new double[grid->row * temp_cols];
			
		}
		else {
			int temp_cols = max<int>(col_process_chunk[col_index + 1], col_process_chunk[col_index]) + 2;
			recv_right_grid = new double[grid->row * temp_cols];
			recv_right_oldgrid = new double[grid->row * temp_cols];
		}

	
		MPI_Irecv(recv_top_grid, 1, grid_bottom, idOfup, 1, MPI_COMM_WORLD, &request[8]);
		MPI_Irecv(recv_bottom_grid, 1, grid_top, idOfdown, 0, MPI_COMM_WORLD, &request[9]);
		MPI_Irecv(recv_left_grid, 1, grid_right, idOfleft, 3, MPI_COMM_WORLD, &request[10]);
		MPI_Irecv(recv_right_grid, 1, grid_left, idOfright, 2, MPI_COMM_WORLD, &request[11]);


		MPI_Irecv(recv_top_oldgrid, 1, grid_bottom, idOfup, 5, MPI_COMM_WORLD, &request[12]);
		MPI_Irecv(recv_bottom_oldgrid, 1, grid_top, idOfdown, 4, MPI_COMM_WORLD, &request[13]);
		MPI_Irecv(recv_left_oldgrid, 1, grid_right, idOfleft, 7, MPI_COMM_WORLD, &request[14]);
		MPI_Irecv(recv_right_oldgrid, 1, grid_left, idOfright, 6, MPI_COMM_WORLD, &request[15]);

		MPI_Waitall(16, request, MPI_STATUS_IGNORE);


		// Assign the row 1 values of recv_ bottom_grid and recv_bottom_oldgrid to the last row of 
		// the grid and old_grid
		if (this->periodic) {
			for (int j = 0; j < grid->col; j++) {
				grid->values[(grid->row - 1) * grid->col + j] = recv_bottom_grid[1 * grid->col + j];
				old_grid->values[(old_grid->row - 1) * old_grid->col + j] = recv_bottom_oldgrid[1 * grid->col + j];
			}
		}
		else {
			// if not periodic and processes with row_index eqaulling p_rows - 1, the recv_bottom_grid 
			// and recv_bottom_oldstri won't receive the values from other process.
			if (row_index != p_rows - 1) {
				for (int j = 0; j < grid->col; j++) {
					grid->values[(grid->row - 1) * grid->col + j] = recv_bottom_grid[1 * grid->col + j];
					old_grid->values[(old_grid->row - 1) * old_grid->col + j] = recv_bottom_oldgrid[1 * grid->col + j];
				}
			}
		}

		// Assign the value of the second-to-last row of recv_top_grid and recv_top_oldgrid to the row 0 
		// of grid  and old_grid
		if (this->periodic) {
			int temp_rows = row_process_chunk[row_index] + 2;
			for (int j = 0; j < grid->col; j++) {
				grid->values[j] = recv_top_grid[(temp_rows - 2) * grid->col + j];
				old_grid->values[j] = recv_top_oldgrid[(temp_rows - 2) * grid->col + j];
			}
		}
		else {
			// if not periodic and processes with row_index eqaulling 0, the recv_top_grid 
			// and recv_top_oldstri won't receive the values from other process.
			if (row_index != 0) {
				int temp_rows = row_process_chunk[row_index] + 2;
				for (int j = 0; j < grid->col; j++) {
					grid->values[j] = recv_top_grid[(temp_rows - 2) * grid->col + j];
					old_grid->values[j] = recv_top_oldgrid[(temp_rows - 2) * grid->col + j];
				}
			}
		}

		// Assign the value of the second-to-last column of recv_left_grid and recv_left_oldgrid to the column 0 
		// of grid and old_grid
		if (this->periodic) {
			int temp_cols = col_process_chunk[col_index] + 2;
			for (int i = 0; i < grid->row; i++) {
				grid->values[i * grid->col] = recv_left_grid[i * temp_cols + col_process_chunk[col_index]];
				old_grid->values[i * grid->col] = recv_left_oldgrid[i * temp_cols + col_process_chunk[col_index]];
			}
		}
		else {
			if (col_index != 0) {
				int temp_cols = col_process_chunk[col_index] + 2;
				for (int i = 0; i < grid->row; i++) {
					grid->values[i * grid->col] = recv_left_grid[i * temp_cols + col_process_chunk[col_index]];
					old_grid->values[i * grid->col] = recv_left_oldgrid[i * temp_cols + col_process_chunk[col_index]];
				}
			}
		}

		// Assign the value of the column 1 of recv_right_grid and recv_right_oldgrid to the last column 
		// of grid and old_grid
		if (this->periodic) {
			int temp_cols = col_process_chunk[col_index] + 2;
			for (int i = 0; i < grid->row; i++) {
				grid->values[i * grid->col+grid->col - 1] = recv_right_grid[i * temp_cols + 1];
				old_grid->values[i * grid->col+ grid->col - 1] = recv_right_oldgrid[i * temp_cols + 1];
			}
		}
		else {
			if (col_index != p_cols - 1) {
				int temp_cols = col_process_chunk[col_index] + 2;
				for (int i = 0; i < grid->row; i ++ ) {
					grid->values[i * grid->col+grid->col-1] = recv_right_grid[i * temp_cols + 1];
					old_grid->values[i * grid->col+grid->col-1] = recv_right_oldgrid[i * temp_cols + 1];
				}
			}
		}
		// do evolution
		this->grid_evolution(new_grid, grid, old_grid, row_index, col_index, p_rows, p_cols,t);
		// write the value to .dat file
		if (t_out <= t) {
			this->grid_to_file(*grid,p_rows, p_cols, row_index, col_index, out_cnt, id);
			out_cnt++;
			t_out += dt_out;
		}
		it++;

		// release memory
		delete[] recv_top_grid; 
		delete[] recv_bottom_grid;
		delete[] recv_left_grid;
		delete[] recv_right_grid;

		delete[] recv_top_oldgrid; 
		delete[] recv_bottom_oldgrid;
		delete[] recv_left_oldgrid;
		delete[] recv_right_oldgrid;

		recv_top_grid = nullptr;
		recv_bottom_grid = nullptr;
		recv_left_grid = nullptr;
		recv_right_grid = nullptr;

		recv_top_oldgrid = nullptr;
		recv_bottom_oldgrid = nullptr;
		recv_left_oldgrid = nullptr;
		recv_right_oldgrid = nullptr;

		MPI_Type_free(&grid_left);
		MPI_Type_free(&grid_right);
		MPI_Type_free(&grid_top);
		MPI_Type_free(&grid_bottom);
		// add Barrier to wait all processes. Then do the next evolution
		MPI_Barrier(MPI_COMM_WORLD);
	}


	if (id == 0) {
		cout << "output: " << out_cnt - 1 << "\tt: " << t << "\titeration: " << it - 1 << endl;
	}

	// release memory
	delete[] request;
	delete new_grid;
	delete grid;
	delete old_grid;
	request = nullptr;
	new_grid = nullptr;
	grid = nullptr;
	old_grid = nullptr;
}


// The evolution of the Grid matrix
void Grid::grid_evolution(Grid*& new_grid, Grid*& grid, Grid*& old_grid, int row_index, int col_index,int p_rows,int p_cols,double& t) {
	for (int i = 1; i < new_grid->row - 1; i++)
		for (int j = 1; j < new_grid->col - 1; j++) {
			new_grid->values[i * grid->col + j] = pow(dt * c, 2.0) * ((grid->values[(i + 1) * grid->col + j] -
				2.0 * grid->values[i * grid->col + j] + grid->values[(i - 1) * grid->col + j]) / pow(dx, 2.0) +
				(grid->values[i * grid->col + j + 1] - 2.0 * grid->values[i * grid->col + j] +
					grid->values[i * grid->col + j - 1]) / pow(dy, 2.0)) + 2.0 * grid->values[i * grid->col + j]
				- old_grid->values[i * grid->col + j];

		}
	
	// If periodic, the boundary value will received from other processes. Therefore, we won't assign values
   
	if(this->Gradient){
		// Assign the left and right columns values
		// left column
		if (col_index == 0) {
			for (int i = 0; i < grid->row; i++) {
				new_grid->values[i * grid->col + 0] = new_grid->values[i * grid->col + 1];
			}
		}

		// right column
		if (col_index == p_cols - 1) {
			for (int i = 0; i < grid->row; i++) {
				new_grid->values[i * grid->col + grid->col - 1] = new_grid->values[i * grid->col + grid->col - 2];
			}
		}

		// assign the top value of the matrix
		if (row_index == 0) {
			for (int j = 0; j < grid->col; j++) {
				new_grid->values[j] = new_grid->values[1 * grid->col + j];
			}
		}

		// assign the bottom value of the matrix
		if (row_index == p_rows - 1) {
			for (int j = 0; j < grid->col; j++) {
				new_grid->values[(grid->row - 1) * grid->col + j] = new_grid->values[(grid->row - 2) * grid->col + j];
			}
		}
	}
	else {
		// The default boundary condition is  Dirichlet
		// the top row  of the matrix
		if (row_index == 0) {
			for (int j = 0; j < grid->col; j++) {
				new_grid->values[j] = 0;
			}
		}
		// the bottom row of the matrix
		if (row_index == p_rows - 1) {
			for (int j = 0; j < grid->col; j++) {
				new_grid->values[(grid->row - 1) * grid->col + j] = 0;
			}
		}

		// the left column
		if (col_index == 0) {
			for (int i = 0; i < grid->row; i++) {
				new_grid->values[i * grid->col] = 0;
			}
		}
		// the right column
		if (col_index == p_cols - 1) {
			for (int i = 0; i < grid->row; i++) {
				new_grid->values[i * grid->col + grid->col - 1] = 0;
			}
		}
		
	}

	t += dt;
	// swap the address to make the grid has the updated value
	Grid* old_address = old_grid;
	old_grid = grid;
	grid = new_grid;
	new_grid = old_address;
}
