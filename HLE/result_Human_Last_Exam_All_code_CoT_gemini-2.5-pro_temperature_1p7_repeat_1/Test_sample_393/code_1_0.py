import math
from collections import Counter

def calculate_tverberg_partitions(a):
    """
    Calculates the number of Tverberg partitions for a set X that always
    induces the numerical partition given by the parts in list 'a'.
    
    The formula is: (r! / (m1! * m2! * ...)) * (n! / (a1! * a2! * ...))
    where:
    - n is the total number of points (sum of a_i)
    - r is the number of parts
    - a_i are the sizes of the parts
    - m_i are the multiplicities of the unique part sizes
    """
    if not a:
        print("The partition list 'a' cannot be empty.")
        return

    # For canonical string representation, sort the input list
    a_sorted = sorted(a)

    # n is the total number of points, r is the number of parts
    n = sum(a_sorted)
    r = len(a_sorted)
    
    print(f"For the given partition of n={n} into r={r} parts: a = {a_sorted}")
    print("-" * 30)
    
    # Calculate factorials of n and r
    n_fact = math.factorial(n)
    r_fact = math.factorial(r)
    
    # --- Denominators ---
    # 1. Product of factorials of the parts a_i
    prod_a_fact = 1
    for part in a_sorted:
        prod_a_fact *= math.factorial(part)
        
    # 2. Product of factorials of the multiplicities of unique part sizes
    counts = Counter(a_sorted)
    multiplicities = sorted(list(counts.values()))
    prod_m_fact = 1
    for m in multiplicities:
        prod_m_fact *= math.factorial(m)
        
    # --- Final Calculation ---
    # We calculate using integer division to avoid float inaccuracies.
    # The formula can be calculated as (r! * n!) // (prod_m_fact * prod_a_fact)
    # The two components are guaranteed to be integers, so we can calculate them separately for clarity.
    
    num_size_permutations = r_fact // prod_m_fact
    num_partitions_for_one_order = n_fact // prod_a_fact
    total_partitions = num_size_permutations * num_partitions_for_one_order

    print("Step 1: Calculate the number of distinct orderings of the part sizes.")
    prod_m_fact_str = "*".join([f"{m}!" for m in multiplicities])
    print(f"   Formula: r! / (m1! * m2! * ...) = {r}! / ({prod_m_fact_str}) = {num_size_permutations}")

    print("\nStep 2: Calculate the number of ways to partition n points for one such ordering.")
    prod_a_fact_str = "*".join([f"{part}!" for part in a_sorted])
    print(f"   Formula: n! / (a1! * a2! * ...) = {n}! / ({prod_a_fact_str}) = {num_partitions_for_one_order}")
    
    print("\nStep 3: The total number is the product of the two results.")
    print("\nThe final equation is:")
    
    # Output the final equation with all numbers
    final_equation = f"({r}! / ({prod_m_fact_str})) * ({n}! / ({prod_a_fact_str})) = {total_partitions}"
    print(final_equation)
    return total_partitions

if __name__ == '__main__':
    # Since no specific partition was given, we use an example.
    # Example partition: n = 7, a_1=3, a_2=2, a_3=2
    example_partition = [3, 2, 2]
    
    calculate_tverberg_partitions(example_partition)
