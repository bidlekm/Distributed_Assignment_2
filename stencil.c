#include "stencil.h"

int main(int argc, char **argv)
{

	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *input = NULL;
	int num_values;
	int num_steps;
	char *output_name;
	char *input_name = argv[1];
	output_name = argv[2];
	num_steps = atoi(argv[3]);
	if (rank == 0)
	{
		if (4 != argc)
		{
			printf("Usage: stencil input_file output_file number_of_applications\n");
			return 1;
		}
		// Read input file
		if (0 > (num_values = read_input(input_name, &input)))
		{
			return 2;
		}
	}

	// Broadcast num_values
	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Scatter input
	int local_num_values = num_values / size;
	double *local_input = malloc(local_num_values * sizeof(double));
	MPI_Scatter(input, local_num_values, MPI_DOUBLE, local_input, local_num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Stencil values
	double h = 2.0 * PI / num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH / 2;
	const double STENCIL[] = {1.0 / (12 * h), -8.0 / (12 * h), 0.0, 8.0 / (12 * h), -1.0 / (12 * h)};
	double leftBoundary[EXTENT];
	double rightBoundary[EXTENT];

	// Allocate data for result
	double *output;
	if (rank == 0)
	{
		if (NULL == (output = malloc(num_values * sizeof(double))))
		{
			perror("Couldn't allocate memory for output in rank 0");
			return 2;
		}
	}

	double *local_output = malloc(local_num_values * sizeof(double));


	MPI_Request leftRequest, rightRequest;
	MPI_Status leftStatus, rightStatus;


	// Start timer
	double start = MPI_Wtime();

	// Repeatedly apply stencil
	for (int s = 0; s < num_steps; s++)
	{
		start -= MPI_Wtime();	
		// Send left and right boundary
		MPI_Isend(&local_input[0], EXTENT, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &leftRequest);
		MPI_Isend(&local_input[local_num_values - EXTENT], EXTENT, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &rightRequest);

		// Receive left and right boundary
		MPI_Irecv(&rightBoundary[0], EXTENT, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, &rightRequest);
		MPI_Irecv(&leftBoundary[0], EXTENT, MPI_DOUBLE, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &leftRequest);
		start += MPI_Wtime();
		for (int i = EXTENT; i < local_num_values - EXTENT; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				int index = i - EXTENT + j;
				result += STENCIL[j] * local_input[index];
			}
			local_output[i] = result;
		}
		start -= MPI_Wtime();
		MPI_Wait(&leftRequest, &leftStatus);
		MPI_Wait(&rightRequest, &rightStatus);
		start += MPI_Wtime();
		for (int i = 0; i < EXTENT; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				int index = i - EXTENT + j;
				double value = (index < 0) ? leftBoundary[EXTENT + index] : local_input[index];
				result += STENCIL[j] * value;
			}
			local_output[i] = result;
		}

		for (int i = local_num_values - EXTENT; i < local_num_values; i++)
		{
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++)
			{
				int index = i - EXTENT + j;
				double value = (index >= local_num_values) ? rightBoundary[index - local_num_values] : local_input[index];
				result += STENCIL[j] * value;
			}
			local_output[i] = result;
		}
		// Swap input and output
		if (s < num_steps - 1)
		{
			double *tmp = local_input;
			local_input = local_output;
			local_output = tmp;
		}
	}

	// Stop timer
	double my_execution_time = MPI_Wtime() - start;

	// Gather output
	MPI_Gather(local_output, local_num_values, MPI_DOUBLE, output, local_num_values, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// // Write result
	// printf("%f\n", my_execution_time);
#ifdef PRODUCE_OUTPUT_FILE
	if (rank == 0 && 0 != write_output(output_name, output, num_values))
	{
		return 2;
	}
#endif

	if (rank == 0)
	{
		printf("%f\n", my_execution_time);
		free(output);
		free(input);
	}


	free(local_output);
	free(local_input);
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