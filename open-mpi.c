#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>

#include <mpi.h>

#define ROOT 0

/*
 * Used for correctness testing
 */
void printMatrix( double *array, int dimen )
{
	int i, j;
	printf( "\n" );
	for( i = 0; i < dimen; i++ )
	{
		for( j = 0; j < dimen; j++ )
		{
			printf( "%.4f\t", array[i*dimen+j] );
		}
		printf( "\n" );
	}
}

/*
 * Load balancing - describes scattering displacement as required by
 * the 'MPI_Scatterv' method to achieve my needs
 */
void spreadDescriber( int *spread_count, int *spread_displs, int p_num,
	int dimen )
{
	int i, index;
	int max_rows, remain_rows;

	int borrow;
	
	if( ( dimen-2 )%p_num != 0 )
	{
		/* this is max num of rows a process will work on */
		max_rows = ( ( dimen-2 ) + p_num - 1 ) / p_num;
	}
	else
	{
		/* or if it divides evenly */
		max_rows = ( dimen-2 )/p_num;
	}
	
	/* total inner rows in the matrix */
	remain_rows = ( dimen-2 );
	
	index = 0; // index to borrow from if needed
	for( i=0; i<p_num; i++ )
	{
		/* workout the the number of */
		if( remain_rows >= max_rows )
		{
			/* + 2*dimen for neighbour rows */
			spread_count[i] = max_rows;
		
			/* decrement the remain rows */
			remain_rows -= max_rows;
		}
		else
		{
			/* if there is no more rows left, then borrow */
			if( remain_rows == 0 )
			{
				/* borrow about half rows from previous p_num */
				borrow = (int) ( spread_count[index] / 2 );
				
				/* spread the numbers back n forth */
				spread_count[i] = borrow;
				spread_count[index] -= borrow;

				/* can't borrow from this index again */
				index++;
			}
			else
			{
				spread_count[i] = remain_rows;
			}
		}
	}
	
	index = 0; // initial displacement postion
	for( i=0; i<p_num; i++ )
	{
		/* times by dimen for row major indexing */
		spread_displs[i] = index;
		index += spread_count[i] *dimen;

		/* now add the actual size and further 2 rows neighbour rows */
		spread_count[i] *= dimen;
		spread_count[i] += 2*dimen;
	}
}

/*
 * Describing the recovery displacement as required by
 * the 'MPI_Gatherv' method to achieve my needs
 */
void gatherDescriber( int *recv_count, int *recv_displs, int *spread_count,
	int *spread_displs, int p_num, int dimen )
{
	int i, items;
	
	/* no. of items to gather from each process */
	for( i=0; i<p_num; i++ )
	{
		items = spread_count[i];

		/* remove the top and bottom rows */
		recv_count[i] = items - 2*dimen;

		/* only interested in the inner row */
		recv_displs[i] = spread_displs[i] + dimen;
	}
}

/*
 * Get the average value from neighbour values
 */
double calculateAverage( double *buff_array, int dimen, int index )
{
	int top, bot, right, left;	//directional index
	int row, col;
	
	row = index / dimen;
	col = index % dimen;
	
	/* row major index workout */
	top = ( row-1 )*dimen + col;
	bot = ( row+1 )*dimen + col;
	right = row*dimen + col+1;
	left = row*dimen + col-1;
	
	return ( 0.25*( buff_array[top] + buff_array[bot] +
			buff_array[right] + buff_array[left] ) );
}

/*
 * Calculate new average values
 */
void updateBuffArrayRows( double *buff_array, double *new_rows, int *spread_count,
	int p_rank, int dimen, int p_num, double *p_max_diff )
{	
	int inner_rows, index;
	double avg_value, current_diff;
	int i,j;
	
	inner_rows = ( spread_count[p_rank]/dimen ) - 2;

	/* creates specific new rows to before replacing the buff array */
	for( i=0; i< inner_rows; i++ )
	{
		index = (i+1)*dimen;	// -> (1...,0) - 2d index
		
		for( j=0; j< dimen; j++ )
		{	
			if( j == 0 )
			{
				/* add current leftmost row value from buff in the row array */
				new_rows[i*dimen] = buff_array[index];
			}
			else if( j == dimen-1 )
			{
				/* add current rightmost row value from buff in the row array */
				new_rows[i*dimen+j] = buff_array[index+j];
			}
			else
			{
				avg_value = calculateAverage(buff_array, dimen, index+j);
				new_rows[i*dimen+j] = avg_value;
				
				current_diff = fabs( avg_value-( buff_array[index+j] ) );
				if( current_diff > *p_max_diff )
				{
					*p_max_diff = current_diff;
				}
			}
		}
	}
	
	/* now update the buff array with the new rows */
	for( i=0; i< inner_rows; i++ )
	{
		index = (i+1)*dimen;
		for( j=0; j< dimen; j++ )
		{
			buff_array[index+j] = new_rows[i*dimen+j];
		}
	}
}

/*
 * Check to see if the dimension and precision is provided.
 * Then see if the provided dimensions and can be divided amongst
 * number of processes
 */
char initialCheck( int p_num, int p_rank, int *dimen, double *precision,
	char *dimen_arg, char *pre_arg )
{
	if( ( !dimen_arg ) | ( !pre_arg ) )
	{
		if( p_rank == ROOT )
		{
			printf( "Dimension or precision not given\n" );
			printf( "E.g. mpirun -np proc_num %s\n",
									"./executable dimension precision" );
		}

		return 'e'; //end
	}

	sscanf( dimen_arg, "%d", dimen );
	sscanf( pre_arg, "%lf", precision );

	if( p_num > ( *dimen-2 ) )
	{
		if( p_rank == ROOT )
		{
			printf( "Too many processes to divide matrix inner rows\n" );
			printf( "given dimen: %d, inner rows: %d, num of p: %d\n",
					*dimen, ( *dimen-2 ), p_num );
		}

		return 'e'; //exit
	}

	return '0';
}

int main(int argc, char **argv)
{
	int p_num, p_rank; 	// no of processes, rank
	
	double *main_matrix;	// significant at root
	double *buff_array;		// buffer chunk of rows from the main array
	
	/* don't really need but helps understanding and debugging */
	int last_inner_row_pos;	// index for first inner row in buff array
	int first_inner_row_pos;// index for last inner row in buff array
	int last_buff_row_pos;	// index for the last row in the buff array
	
	double *new_rows;	// rows to update in buff array
	int new_rows_size;	// its size
	
	int dimen;			// size (dimension) of the array
	double precision;	// precision
	double p_max_diff;	// max difference in current process
	
	double glob_max_diff;	// max difference from all processes
	int iter_count;			// no of iteration

	int i,j; 			// loop counters
	
	/* significant only at the root process to describe scattering */
	int *spread_count;	// no. of items to spread for all the process
	int *spread_displs;	// describing the displacement in scattering the data
	
	/* significant only at the root process to describe recovering */
	int *recv_count;	// no. of items to recover from all the processes
	int *recv_displs;	// describing the displacement to recover the data

	MPI_Status status;
	
	clock_t start, end;	// used for measuring the time
	
	MPI_Init( &argc, &argv );
	
	MPI_Comm_size( MPI_COMM_WORLD, &p_num );
	MPI_Comm_rank( MPI_COMM_WORLD, &p_rank );
	
	// also puts dimension and precision value if check successful
	if( 'e' == initialCheck(
			p_num, p_rank, &dimen, &precision, argv[1], argv[2] ) )
	{
		MPI_Finalize();
		return 0;
	}

	/* array set up by the root process */
	if( p_rank == ROOT )
	{
		/* setupArray( &main_matrix, dimen ); */
		main_matrix = (double *) malloc( dimen*dimen*sizeof( double ) );
		for( i=0; i<dimen; i++ )
		{
			for( j=0; j<dimen; j++ )
			{
				/* make a simple matrix 1 and 0 */
				if( i==0 )
				{
					main_matrix[i*dimen+j] = 1;
				}
				else if( j == 0 )
				{
					main_matrix[i*dimen+j] = 1;
				}
				else
				{
					main_matrix[i*dimen+j] = 0;
				}
				//main_matrix[i*dimen+j] = ( rand() % 9 ) + 1;
			}
		}
		
		//printf( "Initial matrix:\n" );
		//printMatrix( main_matrix, dimen );
	}
	
	/* allocate spread and gather describers based on the no of processes */
	spread_count = malloc( p_num*sizeof( int ) );
	spread_displs = malloc( p_num*sizeof( int ) );
	recv_count = malloc( p_num*sizeof( int ) );
	recv_displs = malloc( p_num*sizeof( int ) );
	
	/* how to spread the array for each processes */
	spreadDescriber( spread_count, spread_displs, p_num, dimen );
	
	/* how and where to put the news from each processes */
	gatherDescriber( recv_count, recv_displs, spread_count, spread_displs,
		p_num, dimen );
	
	buff_array = malloc( spread_count[p_rank]*sizeof( double ) );

	/* allocate space for rows with new values */
	new_rows_size = spread_count[p_rank] - 2*dimen;
	new_rows = malloc( new_rows_size*sizeof( double ) );
	
	if( p_rank == ROOT )
	{
		start = clock();
	}
	
	/* sent out the chunks of rows to each processes */
	MPI_Scatterv( main_matrix, spread_count, spread_displs, MPI_DOUBLE,
					buff_array, spread_count[p_rank], MPI_DOUBLE,
					ROOT, MPI_COMM_WORLD );
	
	/* index in the buff array */
	first_inner_row_pos = dimen;
	last_inner_row_pos = ( spread_count[p_rank]/dimen - 2 )*dimen;
	last_buff_row_pos = last_inner_row_pos + dimen;
	
	iter_count = 0;
	do
	{
		/* reset diff value */
		p_max_diff = 0.0;
		glob_max_diff = 0.0;
		
		/* increment the iteration count */
		iter_count++;
		
		/* get rows with calculated avg values */
		updateBuffArrayRows( buff_array, new_rows, spread_count, p_rank, dimen,
			p_num, &p_max_diff );
		
		/* barrier before process send and receive new rows */
		MPI_Barrier( MPI_COMM_WORLD );
		
		/* send 2nd row to previous process
		 * receive top row from previous process
		 */
		if( p_rank != ROOT )
		{
			MPI_Sendrecv(
				&buff_array[first_inner_row_pos], dimen, MPI_DOUBLE,
					p_rank-1, 1,
				buff_array, dimen, MPI_DOUBLE,
					p_rank-1, 1, MPI_COMM_WORLD, &status );
		}
		
		/* send 2nd last row to next process
		 * receive last row from next process
		 */
		if( p_rank != p_num-1 )
		{
			MPI_Sendrecv(
				&buff_array[last_inner_row_pos], dimen, MPI_DOUBLE,
					p_rank+1, 1,
				&buff_array[last_buff_row_pos], dimen, MPI_DOUBLE,
					p_rank+1, 1,	MPI_COMM_WORLD, &status );	
		}
		
		/* Deduce the max difference held from all the process */
		MPI_Allreduce( &p_max_diff, &glob_max_diff, 1, MPI_DOUBLE, MPI_MAX,
			MPI_COMM_WORLD );
	}
	while( glob_max_diff > precision  ); //loop condition
	
	/* gather the new rows from all the processes back to root */
	MPI_Gatherv( new_rows, recv_count[p_rank], MPI_DOUBLE,
				main_matrix, recv_count, recv_displs, MPI_DOUBLE,
				ROOT, MPI_COMM_WORLD );
	
	if( p_rank == ROOT )
	{	
		end = clock();
		// printf( "\nNew matrix:\n" );
		// printMatrix( main_matrix, dimen );
		
		printf( "Iteration count: %d\n", iter_count );
		printf( "Given precision: %f, %s %f\n", precision,
			"calculation precision at:", glob_max_diff );
		printf( "Dimension: %d, np was: %d", dimen, p_num );
		
		printf( "\nTotal elapsed time: %f\n",
				( double )( end-start )/CLOCKS_PER_SEC );
		
		/* free the allocated main array matrix */
		free( main_matrix );
	}
	
	/* free all the used array */
	free( spread_count );
	free( spread_displs );
	
	free( recv_count );
	free( recv_displs );
	
	free( new_rows );
	free( buff_array );
	
	// MPI finalize
	MPI_Finalize();
	
	return 0;
}