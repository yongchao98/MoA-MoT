import math
from collections import Counter

def solve_tverberg_partition_count():
    """
    This function calculates the number of Tverberg partitions for a given
    set of parameters n and a partition a_1, ..., a_r.

    The problem statement implies that the number of such partitions is a
    combinatorial problem of partitioning a set of n items into unlabeled
    subsets of specified sizes.
    """
    
    # We define a specific instance of the problem to solve.
    # The user can change the list `a` to solve for their own partition.
    # For this example, let's choose a partition a = [4, 3, 3, 2].
    # The geometric constraints (d, etc.) are satisfied by this choice,
    # for instance, with d >= 3.
    a = [4, 3, 3, 2]
    n = sum(a)
    
    print(f"This script calculates the number of Tverberg partitions of a set of n={n} points.")
    print(f"The partition sizes are given by a = {a}.")
    print("-" * 50)
    
    # The formula for the number of ways to partition a set of n distinct items
    # into unlabeled subsets of sizes given by the multiset {a_1, ..., a_r} is:
    # N = n! / ( (a_1! * a_2! * ... * a_r!) * (m_1! * m_2! * ... * m_k!) )
    # where m_j is the number of parts of a certain size.

    # Step 1: Calculate n!
    n_factorial = math.factorial(n)

    # Step 2: Calculate the product of factorials of part sizes (Π a_i!)
    prod_a_factorial = 1
    for part_size in a:
        prod_a_factorial *= math.factorial(part_size)

    # Step 3: Count frequencies of part sizes to find the m_j values
    part_size_counts = Counter(a)
    m_counts = list(part_size_counts.values())
    
    # Step 4: Calculate the product of factorials of these frequencies (Π m_j!)
    prod_m_factorial = 1
    for count in m_counts:
        prod_m_factorial *= math.factorial(count)

    # Step 5: Calculate the final number of partitions
    num_partitions = n_factorial // (prod_a_factorial * prod_m_factorial)
    
    # Step 6: Print the detailed calculation as requested
    print("The number of partitions is calculated using the formula:")
    print("N = n! / ( (Π a_i!) * (Π m_j!) )")
    print("\nStep-by-step calculation:")

    # Print the equation with symbols
    a_fact_str = " * ".join([f"{p}!" for p in sorted(a, reverse=True)])
    m_fact_str = " * ".join([f"{m}!" for m in sorted(m_counts, reverse=True) if m > 1])
    if not m_fact_str: m_fact_str = "1" # Handle case where all parts are unique
    
    print(f"N = {n}! / ( ({a_fact_str}) * ({m_fact_str}) )")
    
    # Print the equation with factorial values
    a_fact_vals_str = " * ".join([str(math.factorial(p)) for p in sorted(a, reverse=True)])
    m_fact_vals_str = " * ".join([str(math.factorial(m)) for m in sorted(m_counts, reverse=True) if m > 1])
    if not m_fact_vals_str: m_fact_vals_str = "1"
    
    print(f"N = {n_factorial} / ( ({a_fact_vals_str}) * ({m_fact_vals_str}) )")
    
    # Print the equation with intermediate products
    print(f"N = {n_factorial} / ( {prod_a_factorial} * {prod_m_factorial} )")
    
    # Print the final fraction
    denominator = prod_a_factorial * prod_m_factorial
    print(f"N = {n_factorial} / {denominator}")
    
    # Print the final result
    print(f"\nThe final number of Tverberg's partitions is: {num_partitions}")
    
    # This return is for capturing the result for the final answer format
    return num_partitions

# Execute the function
result = solve_tverberg_partition_count()
# The final answer is printed above, and we format it below as requested.
# print(f"\n<<< {result} >>>")