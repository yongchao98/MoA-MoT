import math

def calculate_partition_numbers(max_n):
    """
    Calculates the number of partitions p(n) for n from 0 to max_n.
    The number of partitions of n, p(n), is the number of ways of writing n as a sum of positive integers.
    We use a dynamic programming approach.
    p[i] stores the number of partitions of i.
    The update rule is p[i] = p[i] + p[i-k].
    """
    # Initialize a list to store the partition numbers p(n).
    # p[0] is 1 (the empty sum is the only partition of 0).
    partitions = [0] * (max_n + 1)
    partitions[0] = 1

    # Build the table of partition numbers up to max_n.
    for k in range(1, max_n + 1):
        for i in range(k, max_n + 1):
            partitions[i] += partitions[i - k]
            
    return partitions

def solve_problem():
    """
    Solves the problem by calculating the product of d_n for n=1 to 8.
    d_n is the number of partitions of n, p(n).
    """
    N = 8
    
    # Step 1 & 2: Calculate partition numbers p(n) for n up to 8.
    # d_n is the dimension of the fixed subspace, which equals p(n).
    all_partitions = calculate_partition_numbers(N)
    
    # Extract d_n values for n=1 to 8.
    d_values = [all_partitions[n] for n in range(1, N + 1)]

    # Step 3: Compute the product.
    product = 1
    for d_n in d_values:
        product *= d_n
        
    # Step 4: Print the final equation and the result.
    equation_str = " * ".join(map(str, d_values))
    print(f"The values of d_n for n=1 to 8 are: {d_values}")
    print("The product is calculated as:")
    print(f"{equation_str} = {product}")

# Run the solver
solve_problem()
