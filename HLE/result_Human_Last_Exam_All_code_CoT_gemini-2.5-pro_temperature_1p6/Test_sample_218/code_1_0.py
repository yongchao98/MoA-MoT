import math

def partitions_generator(n):
    """
    Generates all integer partitions of n in non-increasing order.
    This is an implementation of Algorithm H from Knuth, TAOCP Vol 4A, Section 7.2.1.4.
    """
    a = [0] * (n + 1)
    k = 1
    a[0] = 0 # sentinal
    a[1] = n
    while k != 0:
        x = a[k] - 1
        y = a[k-1]
        k -= 1
        while x > y:
            a[k+1] = x
            k += 1
            x -= y
        a[k+1] = x + y
        yield a[1:k+2]

def calculate_dimension(n, partition):
    """
    Calculates the dimension of an irreducible representation of S_n
    corresponding to the given partition using the Hook Length Formula.
    """
    # Python's integers handle arbitrary size, so n! is fine for n=25
    n_factorial = math.factorial(n)

    product_of_hooks = 1
    # Iterate over each box (i, j) in the Young diagram
    for i in range(len(partition)): # row index i (0-indexed)
        for j in range(partition[i]): # column index j (0-indexed)
            # Calculate hook length for box (i, j)
            # Number of boxes to the right in the same row
            hooks_right = partition[i] - (j + 1)
            # Number of boxes below in the same column
            hooks_below = 0
            for k in range(i + 1, len(partition)):
                if partition[k] > j:
                    hooks_below += 1
            
            # Total hook length = hooks_right + hooks_below + 1 (for the box itself)
            hook_length = hooks_right + hooks_below + 1
            product_of_hooks *= hook_length

    # Use integer division for precision
    dimension = n_factorial // product_of_hooks
    return dimension

def solve():
    """
    Finds and counts the irreps of S_25 with dimension < 500,000.
    """
    n = 25
    limit = 500000
    count = 0
    
    found_irreps = []

    # Iterate through all partitions of n
    for p in partitions_generator(n):
        # For each partition, calculate the dimension of the corresponding irrep
        dim = calculate_dimension(n, p)
        
        # Check if the dimension is strictly less than the limit
        if dim < limit:
            found_irreps.append((p, dim))
            count += 1
    
    # Sort results by dimension for cleaner output
    found_irreps.sort(key=lambda x: x[1])

    print(f"Found {count} irreducible representations of S_{n} with dimension strictly less than {limit}:")
    print("-" * 80)
    for p, dim in found_irreps:
         print(f"Partition: {str(p):<55} Dimension: {dim}")
    print("-" * 80)
    print(f"Final Count: {count}")

if __name__ == '__main__':
    solve()