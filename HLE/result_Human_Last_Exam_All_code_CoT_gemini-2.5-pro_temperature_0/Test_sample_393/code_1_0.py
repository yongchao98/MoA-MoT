import math
from collections import Counter

def solve_tverberg_partitions():
    """
    Calculates the number of Tverberg partitions for a set X with the properties described.
    The problem is equivalent to counting the number of set partitions of a specific type.
    """
    # This list represents the partition a_1, a_2, ..., a_r.
    # You can change this list to any other valid partition.
    # Example: A partition of n=10 into 4 parts of sizes 3, 3, 2, 2.
    a = [3, 3, 2, 2]

    # Calculate n, the total number of points.
    n = sum(a)

    # Count the occurrences of each part size (for the m_j factorials).
    counts = Counter(a)

    # Calculate the numerator of the formula.
    numerator = math.factorial(n)

    # Calculate the denominator.
    # It is the product of the factorials of the part sizes (a_i!)
    # times the product of the factorials of the counts of identical sizes (m_j!).
    denominator = 1
    denom_parts_str = []

    # Contribution from a_i!
    for part_size in a:
        denominator *= math.factorial(part_size)
        denom_parts_str.append(f"{part_size}!")

    # Contribution from m_j!
    for count in counts.values():
        if count > 1:
            denominator *= math.factorial(count)
            denom_parts_str.append(f"{count}!")

    # The final result is n! / denominator.
    result = numerator // denominator

    # Sort the parts in the denominator string for consistent representation.
    # We sort numerically by extracting the number before the '!'
    denom_parts_str.sort(key=lambda s: int(s[:-1]), reverse=True)

    # Print the final equation as requested.
    print(f"For the partition {a}, the number of Tverberg's partitions is:")
    print(f"{n}! / ({' * '.join(denom_parts_str)}) = {result}")

solve_tverberg_partitions()