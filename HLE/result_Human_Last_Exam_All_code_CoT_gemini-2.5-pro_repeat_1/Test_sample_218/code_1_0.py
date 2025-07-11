import math

def solve():
    """
    Calculates the number of irreducible representations of S_25
    with a dimension strictly less than 500,000.
    """
    n = 25
    threshold = 500000
    count = 0

    # An efficient generator for integer partitions of n.
    # The partition is yielded as a list of integers in descending order.
    def generate_partitions(n):
        a = [0] * (n + 1)
        k = 1
        a[1] = n
        while k != 0:
            x = a[k]
            y = a[k - 1]
            k -= 1
            while x <= y:
                a[k] = x
                y -= x
                k += 1
            a[k] = x + y
            # Yield a copy of the partition slice
            yield a[0:k + 1]

    # Pre-calculate n! as it's used for every dimension calculation.
    try:
        n_factorial = math.factorial(n)
    except OverflowError:
        print(f"Error: {n}! is too large for standard float, but Python handles it.")
        n_factorial = 1
        for i in range(1, n + 1):
            n_factorial *= i


    # Iterate through all partitions of n.
    for p in generate_partitions(n):
        # Calculate the product of hook lengths for the partition p.
        hook_product = 1
        for i in range(len(p)):  # i is the row index (0-based)
            for j in range(p[i]):  # j is the column index (0-based)
                # Calculate the hook length for the cell (i, j).
                
                # Arm: Number of cells to the right in the same row.
                arm = p[i] - (j + 1)
                
                # Leg: Number of cells below in the same column.
                leg = 0
                for k in range(i + 1, len(p)):
                    if p[k] > j:
                        leg += 1
                
                hook_length = arm + leg + 1
                hook_product *= hook_length

        # Calculate the dimension of the irreducible representation.
        # The division n! / hook_product is guaranteed to be an integer.
        dimension = n_factorial // hook_product

        # Check if the dimension is below the threshold.
        if dimension < threshold:
            count += 1
            
    print(f"The number of irreducible representations of S_25 with dimension less than 500,000 is: {count}")

solve()
<<<110>>>