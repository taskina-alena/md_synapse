import numpy as np
import pandas as pd
import sys
from numba import jit


# ##=============================================##
# ##=============================================##
def split_movie(filename, movieFormat, header_lines): ## split dump into time steps

    # if len(header_lines.shape) == 1:
    #     header_lines = header_lines[np.newaxis, :] ## ToDo Single frame read to be imnplemented

    numParts = header_lines[0,2]
    numSplit = len(header_lines)
    start_skip = header_lines[0,1] - 8

    lines_to_skip = []
    for h in header_lines:
        assert h[1] == start_skip+8
        for i in range(start_skip, h[1] + 1):
            lines_to_skip.append(i)
        start_skip = h[1] + numParts + 1


    workingData = np.array(pd.read_csv(filename, skiprows= lines_to_skip, sep = "\s+", names=movieFormat, header = None, dtype=float))

    output = np.reshape( workingData, (numSplit,numParts,workingData.shape[1]))

    assert(output.shape[1] == numParts)

    return output

##=============================================##
##=============================================##
def load_frame(filename, movieFormat, frame):
	read_lines = header_lines(filename)  ## select or slice from read_lines to load single frames ToDo to be implemented

	# print(read_lines[frame])
	# print(read_lines[frame - 1])
	    # start_line = line_instruct[1]
    # max_lines = line_instruct[2]
    # print(line_instruct[1])
    # # out_frame = pd.read_csv('../input/sample_submission.csv', skiprows=5, nrows=10)

    # # out_frame = np.array(pd.read_csv(filename,sep = "\s+", dtype=float, skiprows=start_line, nrows=max_lines))
    # out_frame = np.zeros(shape = (max_lines, len(movieFormat)))
    # counter = 0
    # with open(filename, 'r') as f:
    #     print(counter)
    #     for line in f.readlines()[start_line + 1:start_line + max_lines]:
    #         out_frame[counter] = np.array(line.strip().split(), dtype=float)
    #         counter += 1
    # return out_frame


##=============================================##
##=============================================##
# @jit(nopython = False)
def header_lines(filename):
    """
    Function works only for constant number of partciles. See MC implementation for generalisation.
    :param filename:
    :return:
    """
    counter = 0
    frame_counter = 0
    lines_output = []
    stepper = [0 ,0 ,0]
    with open(filename, 'r') as f: ## pre load first line to determine the size of the output array

        for line in f.readlines():
            line_buffer = line.strip().split()

            if counter == 0: ## some checks on the file format
                if line.startswith("ITEM: TIMESTEP"):
                    pass
                else:
                    print("ERROR: Ooops wrong movie format!")
                    sys.exit()
            if counter == 2:

                if line.startswith("ITEM: NUMBER OF ATOMS"):
                    pass
                else:
                    print("ERROR: Ooops wrong movie format!")
                    sys.exit()
            if counter == 4:

                if line.startswith("ITEM: BOX BOUNDS pp pp pp"):
                    pass
                else:
                    print("ERROR: Ooops wrong movie format!")
                    sys.exit()

            if counter == 3:
                num_parts = int(line_buffer[0])
            if counter == 9:
                num_columns = len(line_buffer)
                # break
            if counter > 0:
                if line.startswith("ITEM: ATOMS"):
                    stepper = [frame_counter, counter , num_parts]
                    lines_output.append(stepper)
                    stepper = [0, 0, 0]
                    frame_counter += 1
            counter += 1
    return np.array(lines_output)
##=============================================##
##=============================================##