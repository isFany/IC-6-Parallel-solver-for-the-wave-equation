import os
import numpy as np
import math
from PIL import Image
import cv2



# The path of the folder
dir_out = os.getcwd() + '\out'
dir_test = os.getcwd() + '\\test'
dir_result = os.getcwd() + '\\result'

# This is the size of the processes. This is same to the C++ file's p
p = 10
# After run the C++ file, the results will show how many files will output.
# Please adjust this parameter same to that value
cnt = 300
# This is the matrix rows. This is same to C++ file's
imax = 100
# This is the matrix columns. This is same to C++ file's
jmax = 100

# This value is used to geneater the vedio. You can change it
fps = 24

# 0 means generate png from the results using stripe decomposition
decomposition = 0


def outPictures(path1, p, cnt, imax, jmax, decomposition):
    """
    Read the data from the .dat file and using the date to output the png file

    Parameters
    ----------
    path1: string
        This is the path of the folder such as os.getcwd()+'\\result'
    p: integer
        The size of the all processes
    cnt: interger
        The number of the .dat file
    imax: integer
        The rows of the matrix
    jmax: integer
        The columns of the matrix
    decomposition: integer
        decide the name of the png. If decoposition is 0, we
        geneater _strip.png, else _grid.png

    Returns
    ----------
           void
    """
    # generate cnt pictures
    for iter in range(0, cnt+1):
        # The .dat file path and the png file path
        if(decomposition == 0):
            arr = np.loadtxt(path1 + '\\result' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_strip.dat')
            outname = 'Pictures\\pic' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_strip.png'
        else:
            arr = np.loadtxt(path1 + '\\result' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_grid.dat')
            outname = 'Pictures\\pic' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_grid.png'
        max_value = arr.max()
        min_value = arr.min()
        array = np.zeros([np.shape(arr)[0], np.shape(arr)[1], 3], dtype=np.uint8)
        # adjust the value in array from 0 to 255
        for i in range(0, np.shape(arr)[0]):
            for j in range(0, np.shape(arr)[1]):
                a=(arr[i][j] - min_value) / (max_value - min_value)*255
                array[i][j]=[0.1*a,0.3*a,0.6*a]
        img = Image.fromarray(array)
        img.save(outname)


def imageTovideo(p, cnt, imax, jmax, fps, decomposition):
    """
    Read the pngs and using these png to output the video

    Parameters
    ----------
    p: integer
        The size of the all processes
    cnt: interger
        The number of the .dat file
    imax: integer
        The rows of the matrix
    jmax: integer
        The columns of the matrix
    fps: integer
        frame rate.
    decomposition: integer
        decide the name of the png. If decoposition is 0, we
        geneater _strip.png, else _grid.png

    Returns
    ----------
           void
    """
    # read the first png file to get the size of these png
    if(decomposition==0):
        img = cv2.imread('Pictures\\pic' + str(0) + '_' + str(imax) + 'x' + str(jmax) + '_strip.png')
    else:
        img = cv2.imread('Pictures\\pic' + str(0) + '_' + str(imax) + 'x' + str(jmax) + '_grid.png')
    imgInfo = img.shape
    size = (imgInfo[1], imgInfo[0])
    videoWriter = cv2.VideoWriter('TestVideo.mp4', -1, fps, size)
    # read cnt number pictures
    for iter in range(0, cnt+1):
        if(decomposition==0):
            frame = cv2.imread('Pictures\\pic' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_strip.png')
        else:
            frame = cv2.imread('Pictures\\pic' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_grid.png')
        videoWriter.write(frame)


def judgeSame(A, B):
    """
    judge whether all values in array A are same to B's or not

    Parameters
    ----------
    A: array
    B: array

    Returns
    ----------
        boolean value. if true, means every value in A is same to B's.
        If false, means some values in A is not same to B's
    """
    for i in range(0, np.shape(A)[0]):
        for j in range(0, np.shape(A)[1]):
            if(abs(A[i][j] - B[i][j]) > 1e-5):
                return False
    return True


def teststrip(path1):
    """
    Compare the results from serial running to the results from parallel
    running to judge whether the results from parallel running is right
    or not. The results are gotten by stripe decomposition

    Parameters
    ----------
    path1: the path of the folder

    Returns
    ----------
        void
    """
    print("----------------Test stripe decomposition start-------------------")
    serial_Gradient1 = np.loadtxt(path1+'\output_Gradient_100x100.dat')
    MPI_Graddient1 = np.loadtxt(path1+'\Grid_MPI_Gradient_100x100.dat')
    result1 = judgeSame(serial_Gradient1, MPI_Graddient1)
    print("Test Neumann boundry size 100x100: ", result1)

    serial_Gradient2 = np.loadtxt(path1 + '\output_Gradient_100x50.dat')
    MPI_Graddient2 = np.loadtxt(path1 + '\Strip_MPI_Gradient_100x50.dat')
    result2 = judgeSame(serial_Gradient2, MPI_Graddient2)
    print("Test Neumann boundry size 100x50: ", result2)

    serial_Dir1 = np.loadtxt(path1 + '\output_Dir_100x100.dat')
    MPI_Dir1 = np.loadtxt(path1 + '\Strip_MPI_Dir_100x100.dat')
    result3 = judgeSame(serial_Dir1, MPI_Dir1)
    print("Test Dirichlet boundry size 100x100: ", result3)

    serial_Dir2 = np.loadtxt(path1 + '\output_Dir_100x50.dat')
    MPI_Dir2 = np.loadtxt(path1 + '\Strip_MPI_Dir_100x50.dat')
    result4 = judgeSame(serial_Dir2, MPI_Dir2)
    print("Test Dirichlet boundry size 100x50: ", result4)

    print("----------------Test stripe decomposition end-------------------")


def testgrid(path1):
    """
    Compare the results from serial running to the results from parallel
    running to judge whether the results from parallel running is right
    or not. The results are gotten by grid decomposition

    Parameters
    ----------
    path1: the path of the folder

    Returns
        void
    """
    print("----------------Test grid decomposition start---------------------")
    serial_Gradient1 = np.loadtxt(path1+'\output_Gradient_100x100.dat')
    MPI_Graddient1 = np.loadtxt(path1+'\Grid_MPI_Gradient_100x100.dat')
    result1 = judgeSame(serial_Gradient1, MPI_Graddient1)
    print("Test Neumann boundry size 100x100: ", result1)

    serial_Gradient2 = np.loadtxt(path1+'\output_Gradient_100x50.dat')
    MPI_Graddient2 = np.loadtxt(path1+'\Grid_MPI_Gradient_100x50.dat')
    result2 = judgeSame(serial_Gradient2, MPI_Graddient2)
    print("Test Neumann boundry size 100x50: ", result2)

    serial_Dir1 = np.loadtxt(path1+'\output_Dir_100x100.dat')
    MPI_Dir1 = np.loadtxt(path1+'\Grid_MPI_Dir_100x100.dat')
    result3 = judgeSame(serial_Dir1, MPI_Dir1)
    print("Test Dirichlet boundry size 100x100: ", result3)

    serial_Dir2 = np.loadtxt(path1+'\output_Dir_100x50.dat')
    MPI_Dir2 = np.loadtxt(path1+'\Grid_MPI_Dir_100x50.dat')
    result4 = judgeSame(serial_Dir2, MPI_Dir2)
    print("Test Dirichlet boundry size 100x100: ", result4)

    print("--------------Test grid decomposition end------------------")


def find_dimension(p):
    """
    Divide the total process p into p_rows × p_cols

    Parameters
    ----------
    p: integer
       The size of the all processes
    
    Returns
    ----------
    p_rows: integer
    p_cols: intger

    Example
    ----------
    >>> m,n = find_dimension(6)
    >>>print(m)
    >>> 2
    >>>print(n)
    >>> 3
    """
    min_gap = p
    top = math.floor(math.sqrt(p) + 1)
    for i in range(1, top + 1):
        if(p % i == 0):
            gap = math.floor(abs(p/i-i))
            if(gap < min_gap):
                min_gap = gap
                p_rows = i
                p_cols = math.floor(p/i)

    return p_rows, p_cols


def setup_partition_grid(p, imax, jmax):
    """
    This function is used to do grid decomposition

    Parameters
    ----------
    p: integer
        The size of the all processes
    imax: integer
        The rows of the matrix
    jmax: integer
        The columns of the matrix


    Returns
    ----------
    row_start: list
            records the row start number in whole matrix for every process
    row_process_chunk: list
            records the number of rows that every process needs to deal with
    col_start: list
           records the column start number in whole matrix for every process
    col_process_chunk: list
           records the the number of colunmns that every process needs to
           deal with
    """
    p_rows, p_cols = find_dimension(p)
    row_start = list(range(p_rows))
    col_start = list(range(p_cols))
    row_process_chunk = list(range(p_rows))
    col_process_chunk = list(range(p_cols))
    # Becasue the first row and last row are boundray, the number of
    # rows to deal with for all processes are the total rows minus 2
    rows_left = imax - 2
    row_start[0] = 0
    for n in range(0, p_rows - 1):
        row_assigned = math.floor(rows_left / (p_rows - n))
        rows_left -= row_assigned
        row_start[n + 1] = row_start[n] + row_assigned
        row_process_chunk[n] = row_assigned
    row_process_chunk[p_rows - 1] = imax - 2 - row_start[p_rows - 1]

    # add the top boundray, process 0's rows start number is still 0
    # the others must add 1
    for i in range(1, p_rows):
        row_start[i] = row_start[i] + 1

    # because the top boundray and the boundary
    # add 1 to the process_chunk
    row_process_chunk[0] = row_process_chunk[0] + 1
    row_process_chunk[p_rows - 1] = row_process_chunk[p_rows - 1] + 1

    # This is same to rows
    col_left = jmax - 2
    col_start[0] = 0
    for n in range(0, p_cols - 1):
        col_assigned = math.floor(col_left / (p_cols - n))
        col_left -= col_assigned
        col_start[n + 1] = col_start[n] + col_assigned
        col_process_chunk[n] = col_assigned
    col_process_chunk[p_cols - 1] = jmax - 2 - col_start[p_cols - 1]

    # add the left boundray and right boundray
    # the others must add 1
    for i in range(1, p_cols):
        col_start[i] = col_start[i] + 1

    # because the left boundary and the last row has the boundary
    # add 1 to the process_chunk
    col_process_chunk[0] = col_process_chunk[0] + 1
    col_process_chunk[p_cols - 1] = col_process_chunk[p_cols - 1] + 1

    return row_start, row_process_chunk, col_start, col_process_chunk


def collectstrip(path1, p, cnt, imax, jmax):
    """
    collect the results from each process using stripe decomposition and
    combine these results to one .dat file

    Parameters
    ----------
    path1:   string
         the path of the folder
    p:  integer
        The size of the all processes
    cnt:  integer
        The number of the .dat file
    imax: integer
        The rows of the matrix
    jamx: integer
        The columns of matrix
    
    Returns
    ----------
        void
    """
    # generate cnt .dat file
    for iter in range(0, cnt+1):
        # the output file path
        outname = 'result\\result' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_strip.dat'
        file = open(outname, 'w')
        for i in range(0, p):
            filename = '\output_' + str(iter) + 'th_' + str(i) + 'of' + str(p) + '_' + str(imax) + 'x' + str(jmax) + '.dat'
            filepath = path1+filename
            for line in open(filepath):
                file.writelines(line)
        file.close()


def collectgrid(path1, p, cnt, imax, jmax):
    """
    collect the results from each process using grid decomposition and
    combine these results to one .dat file

    Parameters
    ----------
    path1:   string
         the path of the folder
    p:  integer
        The size of the all processes
    cnt:  integer
        The number of the .dat file
    imax: integer
        The rows of the matrix
    jamx: integer
        The columns of matrix
    
    Returns
    ----------
        void
    """
    p_rows, p_cols = find_dimension(p)
    # if p_rows, this means teh solver will use stripe decomposition
    # we will use collectstrip too collect the files
    if(p_rows == 1):
        collectstrip(path1, p, cnt, imax, jmax)
    else:
        row_start, row_process_chunk, col_start, col_process_chunk = setup_partition_grid(p, imax, jmax)
        for iter in range(0, cnt+1):
            # initial the array
            Result = np.zeros((imax, jmax))
            outname = 'result\\result' + str(iter) + '_' + str(imax) + 'x' + str(jmax) + '_grid.dat'
            for i in range(0, p):
                row_index = math.floor(i / p_cols)
                col_index = i % p_cols
                fromname ='\output_' + str(iter) + 'th_' + str(i) + 'of('+ str(row_index) + ',' + str(col_index)+')_' + str(imax) + 'x' + str(jmax) + '.dat'
                grid_data = np.loadtxt(path1 + fromname)
                # according the row start number, columns start number,
                # row_process chunk and column process chunk of each process
                # to assign the sub-result to the total result
                for m in range(row_start[row_index], row_start[row_index] + row_process_chunk[row_index]):
                    for n in range(col_start[col_index], col_start[col_index] + col_process_chunk[col_index]):
                        Result[m][n] = grid_data[m - row_start[row_index]][n-col_start[col_index]]
            np.savetxt(outname, Result)


collectgrid(dir_out, p, cnt, imax, jmax)
collectstrip(dir_out, p, cnt, imax, jmax)
teststrip(dir_test)
testgrid(dir_test)
outPictures(dir_result, p, cnt, imax, jmax, decomposition)
imageTovideo(p, cnt, imax, jmax, fps, decomposition)
