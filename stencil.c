#include "stencil.h"

int main(int argc, char **argv)
{
	int rank, size;
    MPI_Status status;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (4 != argc)
	{
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}
	char *input_name = argv[1];
	char *output_name = argv[2];
	int num_steps = atoi(argv[3]);

	// Read input file
	double *input;
	int num_values;

	//TODO: start timer before or after the inputs are read in
	// Start timer
	double start = MPI_Wtime();
	if(rank == 0)
	{
		if (0 > (num_values = read_input(input_name, &input)))
		{
			perror("Couldn't allocate memory for input in process 0");
			return 2;
			//TODO: MPI_Finalize(); and cancel all other threads
		}
		
	}
	// Stencil values
	double h = 2.0 * PI / num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH / 2;
	const double STENCIL[] = {1.0 / (12 * h), -8.0 / (12 * h), 0.0, 8.0 / (12 * h), -1.0 / (12 * h)};
	const int dataPerProcess = num_values / size;
	const int localBufferSize = dataPerProcess + 2*EXTENT;
	double* localBuffer = (double*) malloc(localBufferSize * sizeof(double));
	if(rank == 0)
	{

		//"inner" processes
		for(int i = 1; i < size - 1; ++i)
		{
			MPI_Send(&input[i*dataPerProcess-EXTENT],localBufferSize, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
		}

		//rightmost process buffer
		for(int i = 0; i < localBufferSize; ++i)
		{
			//reuse of buffer: use own local buffer to send data, rewrite with own data later
			localBuffer[i] = input[(((size - 1) * dataPerProcess) + i - EXTENT) % num_values];
			MPI_Send(&localBuffer,localBufferSize, MPI_DOUBLE, size - 1 , 0, MPI_COMM_WORLD);
		}
		//leftmost process buffer (rank 0)
		for(int i = 0; i < localBufferSize; ++i)
		{
			localBuffer[i] = input[(i - EXTENT) % num_values];
		}
	}
	else
	{
    	MPI_Recv(localBuffer, localBufferSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	}
	
	//set up persistent communication with neighboring processes
	MPI_Request receiveRight;
	MPI_Request receiveLeft;
	MPI_Request sendRight;
	MPI_Request sendLeft;
	
	double* sendBuffer = (double*)malloc(EXTENT * sizeof(double));
	double* receiveBuffer = (double*)malloc(EXTENT * sizeof(double));

	MPI_Recv_init(receiveBuffer, EXTENT, MPI_DOUBLE, (rank + size - 1) % size, 0, MPI_COMM_WORLD, &receiveRight);
	MPI_Recv_init(&receiveBuffer[EXTENT], EXTENT, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &receiveLeft);

	MPI_Send_init(sendBuffer, EXTENT, MPI_DOUBLE, (rank + size - 1) % size, 0, MPI_COMM_WORLD, &sendRight);
	MPI_Send_init(&sendBuffer[EXTENT], EXTENT, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &sendLeft);

	// Repeatedly apply stencil
	for (int s = 0; s < num_steps; s++)
	{
		// Apply stencil
		for (int i = 0; i < num_values; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				int index = (i - EXTENT + j + num_values) % num_values;
				result += STENCIL[j] * input[index];
			}
			output[i] = result;
		}
		// Swap input and output
		if (s < num_steps - 1)
		{
			double *tmp = input;
			input = output;
			output = tmp;
		}
	}
	if(rank == 0)
	{
			free(input);
		// Stop timer
		double my_execution_time = MPI_Wtime() - start;

		// Write result
		printf("%f\n", my_execution_time);
#ifdef PRODUCE_OUTPUT_FILE
		if (0 != write_output(output_name, input, num_values))
		{
			return 2;
		}
#endif
	}
	MPI_Finalize();
	return 0;
}

int read_input(const char *file_name, double **values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "r")))
	{
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values))
	{
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double))))
	{
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (EOF == fscanf(file, "%lf", &((*values)[i])))
		{
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int write_output(char *file_name, const double *output, int num_values)
{
	FILE *file;
	if (NULL == (file = fopen(file_name, "w")))
	{
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++)
	{
		if (0 > fprintf(file, "%.4f ", output[i]))
		{
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n"))
	{
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file))
	{
		perror("Warning: couldn't close output file");
	}
	return 0;
}
