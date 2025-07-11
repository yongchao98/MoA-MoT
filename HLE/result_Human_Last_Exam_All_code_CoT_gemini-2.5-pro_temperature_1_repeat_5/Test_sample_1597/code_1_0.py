import math

def solve():
    """
    This function counts the number of "allowed" ordered pairs (a,b) with 1 <= a, b <= 1000.
    A pair is allowed if a=1, or b=1, or if a and b are powers of the same integer g >= 2.
    """
    limit = 1000
    memo_base = {}

    def get_int_root(n, k):
        """Calculates the integer k-th root of n, or returns None if it's not an integer."""
        if k == 1:
            return n
        # Use binary search for a robust integer root finding
        low = 2
        # The root can't be larger than the square root for k>2, or n^(1/k) in general.
        high = int(n**(1/k)) + 1
        while low <= high:
            mid = (low + high) // 2
            if mid == 0:
                low = 1
                continue
            try:
                val = mid**k
            except OverflowError: # val can be very large
                val = float('inf')

            if val == n:
                return mid
            elif val < n and val > 0:
                low = mid + 1
            else:
                high = mid - 1
        return None

    def get_base(n):
        """Finds the 'power base' of n, which is a non-perfect-power g such that n = g^k."""
        if n in memo_base:
            return memo_base[n]

        # The maximum possible exponent k is log2(n)
        max_k = int(math.log2(n))
        for k in range(max_k, 1, -1):
            g = get_int_root(n, k)
            if g is not None:
                # n = g^k. The base of n is the base of g.
                base = get_base(g)
                memo_base[n] = base
                return base

        # If no root is found, n is not a perfect power, so it's its own base.
        memo_base[n] = n
        return n

    # Case 1: Count pairs where a=1 or b=1
    count_ones = limit + limit - 1

    # Case 2: Count pairs where a,b > 1.
    # Group numbers by their power base.
    base_counts = {}
    for n in range(2, limit + 1):
        base = get_base(n)
        base_counts[base] = base_counts.get(base, 0) + 1

    # Calculate the sum of squares of group sizes.
    count_gt_1 = 0
    groups_gt_1_counts = []
    groups_eq_1_count = 0
    for base in sorted(base_counts.keys()):
        count = base_counts[base]
        if count > 1:
            groups_gt_1_counts.append(count)
        else:
            groups_eq_1_count += 1

    for count in groups_gt_1_counts:
        count_gt_1 += count * count
    count_gt_1 += groups_eq_1_count # Add the k=1 cases, where k^2=1

    total_count = count_ones + count_gt_1

    # Print the breakdown of the final calculation
    print(f"The total number of allowed pairs is the sum of two cases:")
    print(f"1. Pairs where a=1 or b=1: {limit} + {limit} - 1 = {count_ones}")
    print(f"2. Pairs where a,b > 1 (sharing a common power base).")
    
    sum_of_squares_parts = [f"{c}^2" for c in groups_gt_1_counts]
    sum_of_squares_parts.append(f"{groups_eq_1_count}*1^2")
    sum_of_squares_str = " + ".join(sum_of_squares_parts)
    
    print(f"   The count for this case is the sum of squares of group sizes: {sum_of_squares_str} = {count_gt_1}")
    print("\nFinal calculation:")
    print(f"{count_ones} + {count_gt_1} = {total_count}")

    print(f"\n<<<{total_count}>>>")

solve()