import math
from collections import Counter

def solve_tverberg_partitions():
    """
    Calculates the number of ways to partition a set of n points into subsets
    of specified sizes, corresponding to the number of Tverberg partitions
    for a specially constructed set X.

    The problem states that for any partition a_1, ..., a_r of n (with 1 <= a_i <= d+1),
    there exists a set X where all Tverberg partitions have this structure.
    This function computes the number of such partitions.
    """

    # --- User-definable input ---
    # Define the partition of n as a list of part sizes a_i.
    # The example below is for a partition of n=11 into parts of size 4, 4, 2, and 1.
    # You can change this list to any other valid partition.
    # For this partition, d must be at least 3 (since max(a) = 4 <= d+1).
    a = [4, 4, 2, 1]

    # --- Calculation ---

    # n is the total number of points, which is the sum of the parts.
    n = sum(a)
    # r is the number of subsets in the partition.
    r = len(a)

    print(f"Calculating the number of partitions for n = {n} with parts a = {a}")
    print("-" * 30)

    # The formula for the number of partitions of a set of n items into
    # unlabeled groups of sizes a_1, ..., a_r is:
    # N = (n!) / (a_1! * ... * a_r!) / (m_1! * ... * m_k!)
    # where m_j is the multiplicity of the j-th distinct part size.

    # 1. Calculate the numerator: n!
    numerator = math.factorial(n)
    print(f"Numerator (n!): {n}! = {numerator}")

    # 2. Calculate the first part of the denominator: a_1! * a_2! * ... * a_r!
    denominator_parts_val = 1
    denominator_parts_str = []
    for part_size in a:
        fact = math.factorial(part_size)
        denominator_parts_val *= fact
        denominator_parts_str.append(f"{part_size}!={fact}")
    print(f"Denominator from parts (Π a_i!): {' * '.join(denominator_parts_str)} = {denominator_parts_val}")


    # 3. Calculate the second part of the denominator for unlabeled groups.
    # This corrects for overcounting when some part sizes are identical.
    # We count the occurrences of each part size to find the multiplicities m_j.
    size_counts = Counter(a) # e.g., {4: 2, 2: 1, 1: 1} for a = [4, 4, 2, 1]

    denominator_multiplicities_val = 1
    denominator_multiplicities_str = []
    for size, count in size_counts.items():
        if count > 1:
            fact = math.factorial(count)
            denominator_multiplicities_val *= fact
            denominator_multiplicities_str.append(f"{count}!={fact}")

    if denominator_multiplicities_val > 1:
        print(f"Part size multiplicities: {dict(size_counts)}")
        print(f"Denominator from multiplicities (Π m_j!): {' * '.join(denominator_multiplicities_str)} = {denominator_multiplicities_val}")
    else:
        print("All part sizes are unique, no multiplicity correction needed.")

    # 4. Calculate the final result.
    total_denominator = denominator_parts_val * denominator_multiplicities_val
    final_result = numerator // total_denominator

    print("-" * 30)
    print("Final Equation:")
    print(f"Number of partitions = {n}! / (({' * '.join([str(p) + '!' for p in a])}) * ({' * '.join([str(c) + '!' for c in size_counts.values() if c > 1]) or '1'}))")
    print(f"                       = {numerator} / ({denominator_parts_val} * {denominator_multiplicities_val})")
    print(f"                       = {numerator} / {total_denominator}")
    print(f"                       = {final_result}")

    return final_result

if __name__ == '__main__':
    number_of_partitions = solve_tverberg_partitions()
    print("\nFinal Answer:")
    print(f"<<<{number_of_partitions}>>>")
