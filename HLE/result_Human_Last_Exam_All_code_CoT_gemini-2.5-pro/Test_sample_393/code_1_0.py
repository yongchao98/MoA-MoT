import math
from collections import Counter

def solve_tverberg_partitions():
    """
    Calculates the number of ways to partition a set of n points into subsets
    with specified sizes, based on the problem description.
    """
    # The problem is posed for "any partition a_1, ..., a_r".
    # As a specific partition is not provided, this script uses an example.
    # You can change the list `a` to any other partition.
    # This example partition a = [2, 2, 1] is valid for Tverberg's theorem
    # with d=1, r=3, n=5.
    a = [2, 2, 1]

    # Calculate the total number of points, n.
    n = sum(a)
    # The number of parts in the partition, r.
    r = len(a)

    print(f"The calculation is for a set of n={n} points with a Tverberg partition structure given by the parts: {a}")
    print("This means we are counting the number of ways to partition a set of {n} items into unlabeled subsets of sizes {a}.")
    print("-" * 30)
    print("The formula is: n! / ( (a_1! * a_2! * ...) * (m_1! * m_2! * ...) )")
    print("where a_i are the part sizes and m_j is the frequency of each part size.")
    print("-" * 30)

    # Numerator: n!
    numerator = math.factorial(n)
    print(f"Step 1: Calculate the numerator n! = {n}! = {numerator}")

    # Denominator Part 1: Product of factorials of part sizes (a_i!)
    denom1 = 1
    denom1_str_list = []
    denom1_val_list = []
    for part_size in a:
        fact_val = math.factorial(part_size)
        denom1 *= fact_val
        denom1_str_list.append(f"{part_size}!")
        denom1_val_list.append(str(fact_val))

    print(f"Step 2: Calculate the product of the factorials of the part sizes:")
    print(f"   {' * '.join(denom1_str_list)} = {' * '.join(denom1_val_list)} = {denom1}")

    # Denominator Part 2: Product of factorials of frequencies (m_s!)
    counts = Counter(a)
    denom2 = 1
    denom2_str_list = []
    denom2_val_list = []
    print(f"Step 3: Count the frequencies of each part size: {dict(counts)}")
    for count in counts.values():
        fact_val = math.factorial(count)
        denom2 *= fact_val
        denom2_str_list.append(f"{count}!")
        denom2_val_list.append(str(fact_val))

    print(f"Step 4: Calculate the product of the factorials of these frequencies:")
    print(f"   {' * '.join(denom2_str_list)} = {' * '.join(denom2_val_list)} = {denom2}")

    # Final Result
    total_denominator = denom1 * denom2
    result = numerator // total_denominator

    print("-" * 30)
    print("Step 5: Combine to find the total number of partitions.")
    final_equation = f"   Number of partitions = {n}! / ( ({' * '.join(denom1_str_list)}) * ({' * '.join(denom2_str_list)}) )"
    print(final_equation)
    print(f"                        = {numerator} / ( {denom1} * {denom2} )")
    print(f"                        = {numerator} / {total_denominator}")
    print(f"                        = {result}")


solve_tverberg_partitions()
<<<15>>>