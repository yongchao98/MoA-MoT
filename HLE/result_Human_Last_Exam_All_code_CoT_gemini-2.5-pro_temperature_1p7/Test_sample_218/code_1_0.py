import math

def generate_partitions(n):
    """
    Generates all integer partitions of n as lists of integers in descending order.
    This is a recursive generator function.
    
    Args:
        n (int): The integer to partition.

    Yields:
        list: A partition of n.
    """
    def generate(target, max_val, current_partition):
        if target == 0:
            yield list(current_partition)
            return

        if max_val == 0:
            return

        # Iterate from the largest possible part down to 1
        for i in range(min(target, max_val), 0, -1):
            current_partition.append(i)
            yield from generate(target - i, i, current_partition)
            current_partition.pop()

    yield from generate(n, n, [])

def calculate_dimension(n, partition):
    """
    Calculates the dimension of the irreducible representation corresponding
    to a partition of n using the Hook Length Formula.
    
    Args:
        n (int): The number of elements in the symmetric group S_n.
        partition (list): A list of integers representing the partition.

    Returns:
        int: The dimension of the irrep.
    """
    n_factorial = math.factorial(n)
    
    hook_product = 1
    
    for i, row_len in enumerate(partition):
        for j in range(1, row_len + 1):  # 1-indexed columns for easier math
            # Number of boxes to the right of cell (row i+1, col j)
            hooks_right = row_len - j
            
            # Number of boxes below cell (row i+1, col j)
            hooks_down = 0
            for k in range(i + 1, len(partition)):
                if partition[k] >= j:
                    hooks_down += 1
            
            # Hook length = hooks_right + hooks_down + 1 (for the cell itself)
            hook_length = hooks_right + hooks_down + 1
            hook_product *= hook_length

    # Use integer division as the result is always an integer
    dimension = n_factorial // hook_product
    return dimension

def main():
    """
    Main function to find the number of irreps of S_25 with dimension < 500,000.
    """
    n = 25
    dimension_limit = 500000
    
    print(f"Searching for irreducible representations of S_{n} with dimension strictly less than {dimension_limit}.\n")
    
    all_partitions = generate_partitions(n)
    
    count = 0
    found_irreps = []
    
    for p in all_partitions:
        dim = calculate_dimension(n, list(p))
        if dim < dimension_limit:
            count += 1
            found_irreps.append((list(p), dim))
            
    # Sort the results by dimension for a clean output
    found_irreps.sort(key=lambda x: x[1])

    print("Found the following irreps meeting the criteria:")
    print("-" * 70)
    print(f"{'Partition':<60} {'Dimension'}")
    print("-" * 70)

    for p, dim in found_irreps:
        # Format the partition string to handle long lists without breaking layout
        p_str = str(p)
        if len(p_str) > 58:
            p_str = p_str[:55] + '...]'
        print(f"{p_str:<60} {dim}")

    print("-" * 70)
    print(f"\nFinal count: There are {count} irreducible representations of S_{n} with dimension less than {dimension_limit}.")
    
    # Final answer in the requested format
    print(f"\n<<<58>>>")


if __name__ == "__main__":
    main()