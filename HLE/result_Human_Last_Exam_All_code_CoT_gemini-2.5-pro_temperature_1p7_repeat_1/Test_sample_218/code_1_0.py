import math

def calculate_dimension(partition):
    """
    Calculates the dimension of an irreducible representation of S_n
    corresponding to a given partition using the Hook Length Formula.
    """
    n = sum(partition)
    
    # Numerator is n!
    try:
        numerator = math.factorial(n)
    except ValueError:
        # n=25 is well within the standard limits, but this handles potential large n.
        return float('inf') 

    # Denominator is the product of all hook lengths
    prod_hook_lengths = 1
    
    # Iterate over each cell (i, j) in the Young diagram of the partition
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            # Calculate the hook length for the cell at (row i, column j)
            # The hook length is 1 (for the cell itself) + arm_length + leg_length.
            
            # arm_length: number of cells to the right in the same row
            arm_length = row_len - (j + 1)
            
            # leg_length: number of cells below in the same column
            leg_length = 0
            for k in range(i + 1, len(partition)):
                # If the row below is long enough to have a cell in this column
                if partition[k] > j:
                    leg_length += 1
            
            hook_length = arm_length + leg_length + 1
            prod_hook_lengths *= hook_length
            
    # The dimension must be an integer, so we use integer division.
    if prod_hook_lengths == 0:
        return float('inf')
    dimension = numerator // prod_hook_lengths
    
    return dimension

def generate_partitions(n, min_part=1):
    """
    A generator for integer partitions of n.
    Yields partitions as lists of integers in non-increasing order.
    """
    if n == 0:
        yield []
        return
        
    for i in range(min_part, n + 1):
        for p in generate_partitions(n - i, i):
            yield [i] + p

def solve_s25_irreps_count():
    """
    Finds and counts the number of irreducible representations of the symmetric group S_25
    that have a dimension strictly less than a specified limit.
    """
    n = 25
    dimension_limit = 500000
    count = 0
    
    # Generate all partitions of 25.
    all_partitions_gen = generate_partitions(n)
    
    # Iterate through all partitions.
    for p in all_partitions_gen:
        # For each partition, calculate the dimension of the corresponding irrep.
        dim = calculate_dimension(p)
        
        # Check if the dimension is within the specified limit.
        if dim < dimension_limit:
            count += 1
    
    print(f"The number of irreducible representations of S_{n} with a dimension strictly less than {dimension_limit} is: {count}")

if __name__ == '__main__':
    solve_s25_irreps_count()
<<<76>>>