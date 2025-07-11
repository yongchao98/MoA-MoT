import math
from collections import Counter

def calculate_number_of_tverberg_partitions(partition_parts):
    """
    Calculates the number of ways to partition a set of n elements into subsets
    with sizes given by the list partition_parts.

    This corresponds to the number of Tverberg partitions for the specific
    point set X described in the problem.

    Args:
        partition_parts (list of int): A list of integers representing the
                                       sizes of the parts in the partition of n.
    """
    if not partition_parts or any(p <= 0 for p in partition_parts):
        print("Error: The partition must contain positive integer parts.")
        return

    # Step 1: Calculate n, the total number of elements.
    n = sum(partition_parts)

    # Step 2: Calculate n!
    n_factorial = math.factorial(n)

    # Step 3: Calculate the product of the factorials of each part size (a_i!).
    prod_a_i_factorial = 1
    for part in partition_parts:
        prod_a_i_factorial *= math.factorial(part)

    # Step 4: Count the multiplicities of each part size (m_s).
    size_counts = Counter(partition_parts)

    # Step 5: Calculate the product of the factorials of these multiplicities (m_s!).
    prod_m_s_factorial = 1
    for count in size_counts.values():
        prod_m_s_factorial *= math.factorial(count)

    # Step 6: Calculate the final number of partitions.
    # N = n! / ( (a_1! * a_2! * ...) * (m_1! * m_2! * ...) )
    if prod_a_i_factorial == 0 or prod_m_s_factorial == 0:
        print("Error: Calculation resulted in a zero denominator.")
        return
        
    number_of_partitions = n_factorial // (prod_a_i_factorial * prod_m_s_factorial)

    # --- Outputting the results clearly ---
    print(f"For the given partition of n = {n} into parts: {partition_parts}")
    print("-" * 30)
    print("The number of Tverberg partitions is the number of ways to partition a set of")
    print(f"{n} elements into groups of sizes {', '.join(map(str, sorted(partition_parts, reverse=True)))}.")

    # Building the string for the final equation
    denom_str_parts = " * ".join(f"{p}!" for p in partition_parts)
    denom_str_mults = " * ".join(f"{c}!" for c, count in size_counts.items() if count > 1)

    print("\nThe formula is: n! / (product of a_i! * product of m_s!)")
    print(f"where n = {n}, a_i are the part sizes, and m_s are their multiplicities.")

    print("\nFinal equation with numbers:")
    
    # We display the multiplicity part only if there are duplicate sizes
    if prod_m_s_factorial > 1:
      counts_of_counts = Counter(size_counts.values())
      mult_fact_str = " * ".join([f"{c}!" for c in size_counts.values() if c > 1])
      
      print(f"Number = {n}! / (({denom_str_parts}) * ({mult_fact_str}))")
      print(f"       = {n_factorial} / (({prod_a_i_factorial}) * ({prod_m_s_factorial}))")
      print(f"       = {number_of_partitions}")
    else: # Simpler case with all unique part sizes
      print(f"Number = {n}! / ({denom_str_parts})")
      print(f"       = {n_factorial} / ({prod_a_i_factorial})")
      print(f"       = {number_of_partitions}")
    print("-" * 30)

if __name__ == '__main__':
    # You can change this list to any partition of n.
    # The problem implies this partition (a_1, ..., a_r) is given.
    # Example: a partition of n=13 into parts of size 4, 4, 2, 2, 1.
    example_partition = [4, 4, 2, 2, 1]
    calculate_number_of_tverberg_partitions(example_partition)
