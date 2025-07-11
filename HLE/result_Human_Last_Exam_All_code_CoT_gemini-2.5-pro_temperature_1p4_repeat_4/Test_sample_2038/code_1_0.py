import math
from collections import Counter
from itertools import permutations

def partitions(n, k, min_val=1):
    """
    Generates all partitions of n into k parts, with each part >= min_val.
    A partition is a tuple of integers in non-decreasing order.
    """
    if k == 1:
        if n >= min_val:
            yield (n,)
        return
    # The first part 'i' can range from min_val up to n/k.
    # The remaining n-i must be partitionable into k-1 parts, each at least i.
    for i in range(min_val, n // k + 1):
        for p in partitions(n - i, k - 1, i):
            yield (i,) + p

def get_knot_count_for_multiset(m):
    """
    Counts the number of unique knots (sequences up to reversal)
    that can be formed from a multiset of parts.
    """
    # Generate all unique permutations of the multiset.
    # Using set to handle duplicates from permutations function for multisets.
    all_perms = set(permutations(m))
    
    palindromes = 0
    non_palindromes_pairs = 0
    
    seen = set()
    
    for p in all_perms:
        if p in seen:
            continue
        
        rev_p = p[::-1]
        
        if p == rev_p:
            palindromes += 1
            seen.add(p)
        else:
            # This pair (p, rev_p) counts as one chiral knot.
            non_palindromes_pairs += 1
            seen.add(p)
            seen.add(rev_p)
            
    return palindromes + non_palindromes_pairs

def solve_and_print():
    """
    Solves the problem and prints the detailed step-by-step calculation.
    """
    total_knots = 0
    max_crossing_number = 13
    
    print("This program calculates the number of 2-bridge knots with crossing number up to 13 that admit two disjoint non-parallel minimal genus Seifert surfaces.")
    print("This is equivalent to counting sequences of positive even integers [c_1, ..., c_n] where n is even and sum(c_i) <= 13, with sequences counted up to reversal.")
    print("-" * 70)
    
    counts_per_C = {}

    # Crossing number C must be a positive even integer.
    for C in range(4, max_crossing_number + 1, 2):
        knots_for_C = 0
        print(f"For crossing number C = {C}:")
        
        # Length of sequence, n, must be even.
        # Smallest part is 2, so max length is C/2.
        for n in range(2, C // 2 + 1, 2):
            # We are partitioning C into n positive even parts.
            # This is equivalent to partitioning C/2 into n parts >= 1.
            half_C = C // 2
            
            # Generate partitions of half_C into n parts.
            for p in partitions(half_C, n, min_val=1):
                multiset = tuple(i * 2 for i in p)
                count = get_knot_count_for_multiset(multiset)
                knots_for_C += count
                
                print(f"  - From partition {multiset} (length {n}): Found {count} knot(s).")

        counts_per_C[C] = knots_for_C
        print(f"Total for C = {C}: {knots_for_C} knot(s)")
        print("-" * 70)
        total_knots += knots_for_C
        
    print("Summary of knots found per crossing number:")
    sum_str_parts = []
    for C in sorted(counts_per_C.keys()):
        print(f"  C={C}: {counts_per_C[C]} knot(s)")
        sum_str_parts.append(str(counts_per_C[C]))

    sum_str = " + ".join(sum_str_parts)
    print(f"\nThe total number of such knots is the sum: {sum_str} = {total_knots}")

if __name__ == '__main__':
    solve_and_print()
    # The final answer in the required format
    print("\n<<<19>>>")
