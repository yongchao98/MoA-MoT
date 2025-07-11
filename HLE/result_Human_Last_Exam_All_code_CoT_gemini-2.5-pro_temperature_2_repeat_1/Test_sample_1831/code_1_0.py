import math
from collections import Counter

def calculate_endomorphism_classes(n):
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size n.
    """
    
    print(f"The task is to find the number of non-isomorphic functional graphs on {n} vertices.")
    print("This is solved by summing counts over all integer partitions of n.\n")

    # u_k is the number of non-isomorphic CONNECTED functional graphs on k vertices.
    # These are pre-computed values from OEIS A000248.
    u_vals = {1: 1, 2: 2, 3: 4, 4: 8, 5: 21}
    if n > 5:
        # These would need to be computed for larger n
        print("This calculation requires pre-computed values for u_k, which are only provided up to k=5.")
        return

    # A simple function to get integer partitions. For n=4, they are easy to list.
    partitions_of_n = {
        1: [[1]],
        2: [[2], [1, 1]],
        3: [[3], [2, 1], [1, 1, 1]],
        4: [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]],
        5: [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]],
    }
    
    partitions = partitions_of_n.get(n)
    if not partitions:
        print(f"Partitions for n={n} are not available.")
        return
        
    def binom_coeff(n_items, k_choices):
        """Calculates the binomial coefficient C(n, k)."""
        return math.comb(n_items, k_choices)

    total_count = 0
    final_equation_parts = []
    print(f"The number of classes is the sum of counts for each partition of {n}:")

    for p in partitions:
        counts = Counter(p)
        term_contribution = 1
        
        # Build text representation for the calculation
        calc_parts = []
        for i in sorted(counts.keys()):
            m_i = counts[i]
            u_i = u_vals[i]
            # Number of ways to choose m_i components of size i from u_i available types
            # is given by the multiset coefficient formula C(u_i + m_i - 1, m_i).
            combinations = binom_coeff(u_i + m_i - 1, m_i)
            term_contribution *= combinations
            calc_parts.append(f"C({u_i}+{m_i}-1, {m_i})")
        
        final_equation_parts.append(str(term_contribution))
        partition_str = " + ".join(map(str, p))
        explanation = f"Partition {partition_str} = {n}:"
        calc_str = " * ".join(calc_parts)
        print(f"- {explanation:<20} ways = {calc_str:<25} = {term_contribution}")
        
        total_count += term_contribution

    print("-" * 60)
    print(f"Total count = {' + '.join(final_equation_parts)} = {total_count}")

# Main execution for the given problem n=4
calculate_endomorphism_classes(4)
<<<18>>>