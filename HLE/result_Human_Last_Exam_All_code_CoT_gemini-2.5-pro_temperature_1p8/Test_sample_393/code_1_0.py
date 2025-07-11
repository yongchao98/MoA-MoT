import math
from collections import Counter

def calculate_tverberg_partitions(a):
    """
    Calculates the number of Tverberg's partitions for a given partition of n.
    This is equivalent to the number of ways to partition a set of n elements
    into unordered subsets of sizes given by the list a.
    
    Args:
        a (list): A list of integers representing the sizes of the parts in the partition.
    """
    
    # Calculate n, the total number of points, by summing the parts.
    n = sum(a)
    
    # Count the frequencies of each part size. This is needed to correct for overcounting
    # when groups of the same size are permuted.
    counts = Counter(a)
    frequencies = list(counts.values())
    
    # Calculate the terms of the formula.
    n_factorial = math.factorial(n)
    
    # Calculate the product of the factorials of the part sizes (denominator of the multinomial coefficient).
    prod_fact_a = 1
    for part in a:
        prod_fact_a *= math.factorial(part)
        
    # Calculate the product of the factorials of the frequencies (the divisor for unordered groups).
    prod_fact_freq = 1
    for freq in frequencies:
        prod_fact_freq *= math.factorial(freq)
        
    # The final result is n! divided by the product of the two denominator parts.
    # We use integer division // to ensure the result is an integer.
    result = n_factorial // (prod_fact_a * prod_fact_freq)
    
    # --- Output the step-by-step calculation ---
    
    print(f"This script calculates the number of ways to partition a set of n points into subsets of given sizes.")
    print(f"The chosen partition of n is: a = {a}")
    print(f"The total number of points is n = sum({a}) = {n}\n")

    print("The formula for the number of partitions of a set of n elements into unordered groups of sizes a_1, a_2, ... is:")
    print("N = (1 / (m_1! * m_2! * ...)) * (n! / (a_1! * a_2! * ...))")
    print("where m_j is the number of times a particular group size appears.\n")

    print("--- Applying the formula to the chosen partition ---\n")
    
    # Display the final equation with numbers filled in.
    final_eq_str_denom_parts = ' * '.join([f"{i}!" for i in a])
    final_eq_str_denom_freqs = ' * '.join([f"{i}!" for i in frequencies])
    
    print(f"Number of Partitions = (1 / ({final_eq_str_denom_freqs})) * ({n}! / ({final_eq_str_denom_parts}))")
    
    # Show the value of each component of the equation.
    print(f"\nBreaking down the calculation:")
    print(f"  n! = {n}! = {n_factorial}")
    
    fact_a_values_str = ' * '.join(map(str, [math.factorial(part) for part in a]))
    print(f"  Product of part factorials ({final_eq_str_denom_parts}) = {fact_a_values_str} = {prod_fact_a}")
    
    freq_values_str = ' * '.join(map(str, [math.factorial(freq) for freq in frequencies]))
    print(f"  Product of frequency factorials ({final_eq_str_denom_freqs}) = {freq_values_str} = {prod_fact_freq}")

    print(f"\nFinal Result:")
    print(f"  Number of Partitions = {n_factorial} / ({prod_fact_a} * {prod_fact_freq})")
    print(f"                       = {n_factorial} / {prod_fact_a * prod_fact_freq}")
    print(f"                       = {result}")


# Since the problem is stated for "any partition", we use a representative
# example to demonstrate the solution. We choose the partition of n=7
# into parts a = [3, 2, 2].
# This corresponds to a setting with d>=2 (since max(a_i) - 1 <= d) and r=3.
example_partition = [3, 2, 2]
calculate_tverberg_partitions(example_partition)
