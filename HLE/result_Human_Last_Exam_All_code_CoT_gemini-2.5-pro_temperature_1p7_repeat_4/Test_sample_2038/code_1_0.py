import math
from itertools import permutations

def partitions(n, k, l=1):
    """
    Generates integer partitions of n into k parts.
    Each part is at least l.
    """
    if k < 1:
        return
    if k == 1:
        if n >= l:
            yield (n,)
        return
    for i in range(l, n // k + 1):
        for part in partitions(n - i, k - 1, i):
            yield (i,) + part

def solve_knot_problem():
    """
    Solves the specific knot theory problem based on its combinatorial equivalent.
    """
    max_crossing_number = 13
    
    # Crossing number S must be an even sum of positive even integers.
    # So we check S = 4, 6, 8, 10, 12.
    
    # S = 2 * B, where B = sum of b_i and a_i = 2 * b_i.
    B_max = max_crossing_number // 2

    knot_counts = {}

    for S in range(4, max_crossing_number + 1, 2):
        B = S // 2
        
        # A set to store unique knot representations for the current crossing number S
        knots_at_S = set()
        
        # n is the number of coefficients, must be even.
        # The smallest parts can be 1, so n can't be larger than B.
        for n in range(2, B + 1, 2):
            # Find partitions of B into n parts. Each part corresponds to a b_i.
            for part in partitions(B, n):
                # Generate unique compositions from this partition.
                for p in set(permutations(part)):
                    # The sequence of even numbers is [2*b_1, 2*b_2, ...]
                    a = tuple(2 * x for x in p)
                    
                    # A knot and its mirror are non-distinct. The knot from sequence `a`
                    # is the mirror of the knot from the reverse of `a`.
                    # We use a frozenset to count the pair {a, reversed(a)} as one item.
                    a_rev = tuple(reversed(a))
                    knot_representation = frozenset([a, a_rev])
                    knots_at_S.add(knot_representation)
        
        count_at_S = len(knots_at_S)
        if count_at_S > 0:
            knot_counts[S] = count_at_S

    print("The number of 2-bridge knots with crossing number at most 13 admitting two disjoint non-parallel minimal genus Seifert surfaces is determined by counting sequences of positive even integers.")
    print("The counts for each crossing number are:")

    counts = []
    for S in sorted(knot_counts.keys()):
        count = knot_counts[S]
        counts.append(count)
        print(f"Crossing number {S}: {count} knot(s)")

    total_knots = sum(counts)
    
    sum_str = " + ".join(map(str, counts))
    
    print(f"\nTotal: {sum_str} = {total_knots}")

if __name__ == "__main__":
    solve_knot_problem()
