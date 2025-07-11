import math
from collections import Counter

def count_set_partitions_with_sizes(a):
    """
    Calculates the number of ways to partition a set of n elements into
    unlabeled subsets of specified sizes.

    The formula is: n! / ( (product of a_i!) * (product of m_j!) )
    where a_i are the part sizes and m_j is the frequency of each part size.
    """
    if not isinstance(a, list) or not all(isinstance(x, int) and x > 0 for x in a):
        print("Error: Input must be a list of positive integers.")
        return

    # Calculate n, the total number of elements
    n = sum(a)

    # Calculate the numerator: n!
    numerator = math.factorial(n)

    # Calculate the first part of the denominator: product of a_i!
    prod_a_i_factorial = 1
    for part_size in a:
        prod_a_i_factorial *= math.factorial(part_size)

    # Calculate the second part of the denominator: product of m_j!
    # where m_j is the number of parts with the same size.
    counts = Counter(a)
    prod_m_j_factorial = 1
    for size in counts:
        m_j = counts[size]
        prod_m_j_factorial *= math.factorial(m_j)
    
    denominator = prod_a_i_factorial * prod_m_j_factorial
    
    # Calculate the final result
    result = numerator // denominator

    # Print the detailed breakdown of the calculation as requested.
    print(f"The calculation is for a partition of a set of n points into parts of sizes a_i.")
    print(f"Given partition sizes: a = {a}")
    print(f"Total number of points: n = sum(a) = {n}")
    
    den_part1_str = " * ".join([f"{x}!" for x in a])
    den_part2_str = " * ".join([f"{counts[size]}!" for size in sorted(counts.keys())])
    
    print("\nThe number of partitions is given by the formula:")
    print(f"N = n! / ( (a_1! * a_2! * ...) * (m_1! * m_2! * ...) )")
    print(f"N = {n}! / ( ({den_part1_str}) * ({den_part2_str}) )")
    
    print("\nPlugging in the numbers:")
    print(f"N = {numerator} / ( ({prod_a_i_factorial}) * ({prod_m_j_factorial}) )")
    print(f"N = {numerator} / {denominator}")
    print(f"N = {result}")

    return result

# Based on the problem description, we can construct a valid example.
# Let's choose dimension d=2 and number of parts r=4.
# Tverberg's theorem applies for n >= (r-1)(d+1)+1 = (4-1)(2+1)+1 = 10.
# The part sizes must be 1 <= a_i <= d+1 = 3.
# A valid partition of n=10 into r=4 parts with a_i <= 3 is 3+3+2+2.
example_partition = [3, 3, 2, 2]

# Calculate and print the number of Tverberg partitions for this example.
final_answer = count_set_partitions_with_sizes(example_partition)
