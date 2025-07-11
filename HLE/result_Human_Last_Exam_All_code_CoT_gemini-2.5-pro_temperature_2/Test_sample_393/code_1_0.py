import math
from collections import Counter

def calculate_tverberg_partitions():
    """
    Calculates the number of Tverberg's partitions for a given partition of n.
    The problem asks for the number of ways to partition a set of n points,
    given that every such partition consists of subsets of specified sizes a_1, ..., a_r.

    This translates to a combinatorial problem: counting the number of ways to partition
    a set of n labeled items into unlabeled subsets of given sizes.

    The formula is:
      Number of Partitions = ( Multinomial Coefficient ) / ( Product of factorials of counts of same-sized parts )
                           = ( n! / (a_1! * a_2! * ...) ) / ( m_1! * m_2! * ... )
    
    This script will calculate this value for a specific example.
    """
    
    # --- Example Input ---
    # We choose a non-trivial example where some parts have the same size
    # to demonstrate the full formula.
    # Let n = 10 be the total number of points.
    # Let the partition of the integer n be given by the parts [3, 3, 2, 2].
    # This means we are partitioning 10 points into four groups: two of size 3, and two of size 2.
    n = 10
    parts = [3, 3, 2, 2]

    # --- Step 1: Calculate the multinomial coefficient part ---
    # Numerator is n!
    numerator_val = math.factorial(n)
    
    # Denominator is the product of the factorials of the part sizes (a_i!)
    denominator_parts_val = 1
    for p in parts:
        denominator_parts_val *= math.factorial(p)
    
    # This gives the number of ways to form *ordered* groups
    multinomial_coeff = numerator_val // denominator_parts_val

    # --- Step 2: Calculate the correction factor for indistinguishable groups ---
    # Count the frequency of each part size (e.g., size 3 appears twice)
    part_size_counts = Counter(parts)
    
    # The divisor is the product of factorials of these frequencies (m_s!)
    # For parts [3, 3, 2, 2], the counts are {3: 2, 2: 2}. Divisor is 2! * 2!.
    indistinguishability_divisor_val = 1
    for count in part_size_counts.values():
        indistinguishability_divisor_val *= math.factorial(count)

    # The final result is the number of partitions into unordered (indistinguishable) groups
    num_partitions = multinomial_coeff // indistinguishability_divisor_val
    
    # --- Step 3: Format the output equation string as requested ---
    # Build the string for the denominator of the multinomial coefficient
    den_parts_str = " * ".join([f"{p}!" for p in parts])
    
    # Build the string for the divisor for indistinguishable groups
    # This only includes counts greater than 1
    count_strs = [f"{c}!" for c in part_size_counts.values() if c > 1]
    den_counts_str = " * ".join(count_strs)
    
    # Combine everything into the final equation string
    equation_str = f"({n}! / ({den_parts_str})) / ({den_counts_str}) = {num_partitions}"

    print(f"The number of Tverberg's partitions for a set of n={n} points inducing the integer partition {parts} is calculated as follows:")
    print(equation_str)


# Execute the calculation
calculate_tverberg_partitions()