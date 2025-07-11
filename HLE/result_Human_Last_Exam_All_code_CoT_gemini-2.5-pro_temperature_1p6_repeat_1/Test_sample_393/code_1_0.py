import math
from collections import Counter

def solve_tverberg_partition_count():
    """
    Calculates the number of Tverberg partitions for a set of points X
    constructed to have all its Tverberg partitions conform to a specific
    size distribution.

    The problem boils down to a combinatorial counting problem: "How many ways
    can a set of n labeled items be partitioned into unlabeled subsets of
    given sizes a_1, a_2, ..., a_r?"

    The formula is: n! / ( (a_1! * ... * a_r!) * (k_1! * ... * k_m!) )
    where:
    - n is the total number of points.
    - a_i are the sizes of the subsets in the partition.
    - k_j is the frequency of the j-th distinct size among the a_i's.

    Since the prompt does not provide a specific partition, we will use a
    concrete example: a partition of n=9 points into 4 subsets of sizes
    [3, 2, 2, 2]. This is a valid partition for d>=2, as max(a_i) <= d+1.
    """

    # Example partition a = [a_1, a_2, ..., a_r].
    partition_sizes = [3, 2, 2, 2]

    # Calculate n, the total number of points.
    n = sum(partition_sizes)

    # Calculate the numerator of the formula: n!
    numerator = math.factorial(n)

    # Calculate the first part of the denominator: product of factorials of part sizes.
    denom_part1 = 1
    for size in partition_sizes:
        denom_part1 *= math.factorial(size)

    # Count the frequencies of each distinct part size.
    # For [3, 2, 2, 2], the counts are {3: 1, 2: 3}.
    size_counts = Counter(partition_sizes)

    # Calculate the second part of the denominator: product of factorials of frequencies.
    # This corrects for overcounting when groups have the same size.
    denom_part2 = 1
    for count in size_counts.values():
        denom_part2 *= math.factorial(count)

    # Calculate the final number of partitions.
    # We use integer division // as the result must be an integer.
    total_denominator = denom_part1 * denom_part2
    num_partitions = numerator // total_denominator
    
    # As requested, output the components of the final equation.
    print(f"The problem is to find the number of ways to partition a set of n points according to the given sizes.")
    print(f"Let's use the example partition a = {partition_sizes}.")
    print("-" * 30)
    
    print(f"1. The total number of points is n = sum({partition_sizes}) = {n}.")
    
    print(f"\n2. The numerator is n! = {n}! = {numerator}.")
    
    denom_part1_str = " * ".join([f"{s}!" for s in partition_sizes])
    print(f"\n3. The first part of the denominator is the product of the factorials of the part sizes:")
    print(f"   {denom_part1_str} = {denom_part1}")

    counts_str = ", ".join([f"size {s} appears {c} time(s)" for s, c in size_counts.items()])
    denom_part2_str = " * ".join([f"{c}!" for c in size_counts.values()])
    print(f"\n4. The frequencies of the part sizes are: {counts_str}.")
    print(f"   To correct for indistinguishable groups of the same size, we use the second part of the denominator:")
    print(f"   {denom_part2_str} = {denom_part2}")
    
    print("-" * 30)
    print(f"The final calculation for the number of partitions is:")
    print(f"   Number = {n}! / (({denom_part1_str}) * ({denom_part2_str}))")
    print(f"          = {numerator} / (({denom_part1}) * ({denom_part2}))")
    print(f"          = {numerator} / {total_denominator}")
    print(f"          = {num_partitions}")

solve_tverberg_partition_count()
