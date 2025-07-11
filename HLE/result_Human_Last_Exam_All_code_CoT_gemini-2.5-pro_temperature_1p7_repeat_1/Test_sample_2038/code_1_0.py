import math
from itertools import permutations

def modInverse(a, m):
    """Calculates the modular multiplicative inverse of a modulo m."""
    if math.gcd(a, m) != 1:
        raise ValueError("Modular inverse does not exist")
    return pow(a, -1, m)

def compute_pq(sequence):
    """
    Computes p, q for the continued fraction p/q = [a_1, a_2, ..., a_n].
    """
    if not sequence:
        return 0, 1
    # Initialize with the fraction a_n / 1
    p, q = sequence[-1], 1
    # Iterate backwards from the second to last element
    for i in range(len(sequence) - 2, -1, -1):
        # p_new/q_new = a_i + 1/(p/q) = a_i + q/p = (a_i*p + q) / p
        p, q = sequence[i] * p + q, p
    return p, q

def find_partitions(n):
    """
    Finds all integer partitions of n. Yields lists in descending order.
    """
    if n == 0:
        yield []
        return
    for i in range(n, 0, -1):
        for p in find_partitions(n - i):
            if not p or i <= p[0]:
                yield [i] + p

def solve_and_print():
    """
    Solves the knot problem and prints the results as per the user's request.
    """
    max_crossing = 13
    knot_details = {} # Store details to avoid recounting: repr -> (c, sequence)

    # Crossing number `c` must be even. Smallest is 4 ([2,2]).
    for c in range(4, max_crossing + 1, 2):
        c_half = c // 2
        
        # Find partitions of c/2. This gives sets of {k_i}.
        # The even coefficients of the continued fraction are {2*k_i}.
        for part in find_partitions(c_half):
            
            # For p to be odd, the length of the sequence `n` must be even.
            if len(part) % 2 != 0:
                continue

            a_sequence_multiset = [i * 2 for i in part]
            
            # Generate unique permutations of this multiset
            for seq_tuple in set(permutations(a_sequence_multiset)):
                sequence = list(seq_tuple)
                
                # Calculate p and q from the continued fraction
                p_val, q_val = compute_pq(sequence)
                
                # Get the canonical representation of the knot K(p,q)
                try:
                    q_inv = modInverse(q_val, p_val)
                    orbit = {q_val, p_val - q_val, q_inv, p_val - q_inv}
                    canonical_q = min(orbit)
                    knot_repr = (p_val, canonical_q)
                    
                    if knot_repr not in knot_details:
                         knot_details[knot_repr] = c
                except ValueError:
                    # Should not happen since gcd(p,q)=1
                    continue
    
    # Process results for printing
    counts = {}
    for c_val in knot_details.values():
        counts[c_val] = counts.get(c_val, 0) + 1
        
    total_knots = len(knot_details)
    
    print("Number of 2-bridge knots with crossing number C <= 13 admitting two disjoint non-parallel minimal genus Seifert surfaces:")
    
    equation_parts = []
    for c in sorted(counts.keys()):
        print(f"- {counts[c]} knot(s) with C = {c}")
        equation_parts.append(str(counts[c]))

    final_equation = " + ".join(equation_parts)
    print("\nFinal Calculation:")
    print(f"Total = {final_equation} = {total_knots}")

    # The final answer in the requested format
    print(f"\n<<<{total_knots}>>>")


solve_and_print()