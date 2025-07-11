import math
from collections import Counter

def get_partitions_recursive(n, min_part=1):
    """
    Generates all integer partitions of n.
    """
    if n == 0:
        yield []
        return
    for i in range(min_part, n + 1):
        for p in get_partitions_recursive(n - i, i):
            yield [i] + p

def get_centralizer_size(partition):
    """
    Calculates the size of the centralizer z(lambda) for a given partition lambda.
    """
    counts = Counter(partition)
    z = 1
    for k, c_k in counts.items():
        z *= (k**c_k) * math.factorial(c_k)
    return z

def get_conjugacy_class_size(partition, n):
    """
    Calculates the size of the conjugacy class for a given partition.
    """
    if sum(partition) != n:
        raise ValueError("Partition must sum to n")
    z = get_centralizer_size(partition)
    return math.factorial(n) // z

def solve():
    """
    Solves the problem for m=3.
    """
    m = 3
    
    # Step 1: Determine n and the partition lambda
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    lambda_partition = []
    for k in range(1, m + 1):
        lambda_partition.extend([k] * (m + 1 - k))
    lambda_partition.sort(reverse=True)

    # Step 2: Calculate the components of the formula
    
    # Component 1: Size of the conjugacy class C_lambda
    size_C_lambda = get_conjugacy_class_size(lambda_partition, n)
    
    # Component 2: n!
    n_factorial = math.factorial(n)
    
    # Component 3: Sum of centralizer sizes over all partitions of n
    all_partitions_of_n = get_partitions_recursive(n)
    sum_of_centralizer_sizes = sum(get_centralizer_size(p) for p in all_partitions_of_n)

    # Step 3: Calculate the final cardinality
    # Formula: |C_lambda| * n! * (sum of z(lambda'))
    cardinality = size_C_lambda * n_factorial * sum_of_centralizer_sizes
    
    # Print the components of the final equation
    print(f"The final equation is |C_lambda| * n! * (sum of z(lambda'))")
    print(f"|C_lambda| = {size_C_lambda}")
    print(f"n! = {n_factorial}")
    print(f"Sum of z(lambda') for all partitions of n={n} is {sum_of_centralizer_sizes}")
    
    # Print the final result
    print("\nFinal Result:")
    print(f"The cardinality of Theta_inverse(lambda) is: {cardinality}")
    
    # Print the first 40 digits of the result
    cardinality_str = str(cardinality)
    first_40_digits = cardinality_str.ljust(40, '0')[:40]
    print(f"\nThe first 40 digits of the cardinality are:")
    print(first_40_digits)

solve()