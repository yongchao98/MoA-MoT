import collections

# Memoization cache for the partition function to improve performance.
partition_memo = {}

def get_integer_partitions(n):
    """
    Generates all unique partitions of an integer n.
    Each partition is a list of integers that sum to n.
    The result is cached using memoization.
    """
    if n in partition_memo:
        return partition_memo[n]
    if n == 0:
        return [[]]
    if n < 0:
        return []

    partitions = []
    # Using a standard recursive algorithm to find partitions.
    # To avoid duplicates, we enforce that parts are in non-increasing order.
    def find_partitions_recursive(target, max_part, current_partition):
        if target == 0:
            partitions.append(current_partition)
            return
        if target < 0:
            return

        for i in range(min(target, max_part), 0, -1):
            find_partitions_recursive(target - i, i, current_partition + [i])

    find_partitions_recursive(n, n, [])
    partition_memo[n] = partitions
    return partitions

def solve_knot_problem():
    """
    Solves the specified knot theory problem by counting integer partitions
    with specific properties.
    """
    print("Finding the number of 2-bridge knots with crossing number at most 13")
    print("that have at least two disjoint non-parallel minimal genus Seifert surfaces.\n")
    print("This property corresponds to knots whose continued fraction has only positive, even coefficients,")
    print("and the number of coefficients is even.\n")

    total_knots = 0
    counts_per_crossing = []

    # The crossing number C is a sum of positive even integers, so C must be even.
    # The smallest crossing number for a knot is 3.
    for c in range(4, 14):
        if c % 2 != 0:
            continue

        # Let C = a_1 + ... + a_k, where a_i are positive and even, and k is even.
        # Let N = C/2. Then N = b_1 + ... + b_k, where b_i = a_i/2 are positive integers.
        # We need to find the number of partitions of N into an even number of parts.
        n = c // 2
        
        all_partitions_of_n = get_integer_partitions(n)
        
        valid_partitions_of_n = []
        for p in all_partitions_of_n:
            if len(p) % 2 == 0:
                valid_partitions_of_n.append(p)
        
        count_for_c = len(valid_partitions_of_n)
        
        if count_for_c > 0:
            print(f"For crossing number C = {c}:")
            print(f"  This corresponds to partitions of N = {c//2} into an even number of parts.")
            
            for p in valid_partitions_of_n:
                # The partition of N gives b_i, we need a_i = 2*b_i for the continued fraction.
                cf_partition = [i * 2 for i in p]
                print(f"  - Partition of {n}: {p}  ->  Continued Fraction terms (partition of {c}): {cf_partition}")
            
            print(f"  Number of knots found for C={c}: {count_for_c}\n")
            counts_per_crossing.append(count_for_c)
            total_knots += count_for_c

    print("--------------------------------------------------")
    print("Summary:")
    print("The number of such knots for each crossing number is:")
    crossing_numbers = [c for c in range(4, 14) if c % 2 == 0]
    for c, count in zip(crossing_numbers, counts_per_crossing):
         print(f"  C = {c}: {count} knot(s)")

    # Final equation as requested.
    equation_parts = [str(x) for x in counts_per_crossing]
    equation_str = " + ".join(equation_parts)
    print(f"\nTotal number of such knots = {equation_str} = {total_knots}")

solve_knot_problem()
<<<14>>>