import math

def solve_power_series_problem():
    """
    This function analyzes the properties of four sets of natural numbers to determine
    for which of them a specific type of power series can exist.

    The problem asks for which subsets S of natural numbers there exists a power series
    with non-zero coefficients only at indices in S, such that the series converges
    everywhere on the closed unit disc (|z| <= 1), but does not converge absolutely
    for |z| = 1.

    This property is directly linked to a concept in harmonic analysis called a "Sidon set".
    A power series with the desired properties exists if and only if the set of indices S
    is NOT a Sidon set.

    We will analyze each set based on known mathematical results.
    """

    # --- Analysis of each set ---

    # Set 1: S = {sum_{k<=n} N_k} where N_k ~ Poi(1) (the range of a 1D random walk)
    # Analysis: The range of a one-dimensional random walk with non-negative steps
    # (like the one described) is known to have significant additive structure.
    # This structure almost surely prevents it from being a Sidon set.
    is_sidon_1 = False  # Almost surely not a Sidon set.

    # Set 2: S = {n^k} for k >= 4
    # Analysis: It is a known result in harmonic analysis that the set of k-th powers
    # of integers, {n^k}, is not a Sidon set for any integer k >= 2. Such sets
    # are "convex" sequences, and this property is related to them not being Sidon sets.
    is_sidon_2 = False # Not a Sidon set.

    # Set 3: The set of primes
    # Analysis: It was a famous open problem for a long time, but it has been proven
    # that the set of prime numbers IS a Sidon set. This is a deep result related
    # to the lack of additive structure in the primes.
    is_sidon_3 = True  # Is a Sidon set.

    # Set 4: S = {floor((pi/2)^n)}
    # Analysis: Let q = pi/2 ≈ 1.57 > 1. The elements s_n = floor(q^n) grow exponentially.
    # The ratio of consecutive terms s_{n+1}/s_n approaches q, which is greater than 1.
    # This means the set is a "Hadamard gap set". All Hadamard gap sets are Sidon sets.
    is_sidon_4 = True  # Is a Sidon set.

    # --- Conclusion ---
    # The question asks for which sets the desired power series can exist.
    # This is true if and only if the set is NOT a Sidon set.

    sets_with_property = []
    if not is_sidon_1:
        sets_with_property.append(1)
    if not is_sidon_2:
        sets_with_property.append(2)
    # The following would be added if they were not Sidon sets
    # if not is_sidon_3:
    #     sets_with_property.append(3)
    # if not is_sidon_4:
    #     sets_with_property.append(4)

    print("A power series with the specified properties exists if and only if the set of indices is NOT a Sidon set.")
    print("-" * 30)
    print(f"Set 1 (Random walk range) is a Sidon set: {is_sidon_1} (almost surely)")
    print(f"Set 2 (Integer powers k>=4) is a Sidon set: {is_sidon_2}")
    print(f"Set 3 (Primes) is a Sidon set: {is_sidon_3}")
    print(f"Set 4 (Exponential sequence) is a Sidon set: {is_sidon_4}")
    print("-" * 30)
    
    # We don't have a numerical equation, but we can print the final logical conclusion.
    # The problem asks for the numbers of the sets that satisfy the condition.
    print(f"The sets that satisfy the condition are: {sets_with_property[0]} and {sets_with_property[1]}")


solve_power_series_problem()