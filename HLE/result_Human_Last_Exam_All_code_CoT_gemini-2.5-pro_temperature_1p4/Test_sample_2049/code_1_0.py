import math
from collections import Counter

def solve():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction outside the prime 2.

    This is equivalent to counting the number of quintic étale algebras over Q
    unramified outside 2.
    """
    # N[d] is the number of number fields of degree d over Q unramified outside 2.
    # These values are obtained from number theory databases (e.g., LMFDB).
    N = {1: 1, 2: 3, 3: 0, 4: 7, 5: 0}

    # All partitions of the integer 5 represent the possible structures of a
    # quintic étale algebra as a product of number fields.
    partitions_of_5 = [
        [5],
        [4, 1],
        [3, 2],
        [3, 1, 1],
        [2, 2, 1],
        [2, 1, 1, 1],
        [1, 1, 1, 1, 1]
    ]

    print("Calculating the number of algebras for each partition of 5:")
    print("-" * 50)

    total_count = 0
    term_counts = []

    # For each partition, we count the number of ways to form an algebra.
    for part in partitions_of_5:
        # Get multiplicities of each degree in the partition
        multiplicities = Counter(part)
        
        term_result = 1
        
        # Calculate the number of ways for this partition structure using the formula:
        # Product over d of C(N_d + m_d - 1, m_d)
        # where m_d is the multiplicity of d in the partition.
        for degree in sorted(multiplicities.keys()):
            m_d = multiplicities[degree]
            N_d = N.get(degree, 0)
            term_result *= math.comb(N_d + m_d - 1, m_d)
        
        # Print the breakdown for this partition
        part_str = " + ".join(map(str, part))
        print(f"Partition {part_str}:")
        
        calc_str_parts = []
        for degree in sorted(multiplicities.keys()):
            m_d = multiplicities[degree]
            N_d = N.get(degree, 0)
            if N_d > 0 or m_d > 1:
                calc_str_parts.append(f"C({N_d}+{m_d}-1, {m_d})")
            else: # C(0,1)
                calc_str_parts.append(f"N_{degree}")
        
        calc_str = " * ".join(calc_str_parts)
        print(f"  Count = {calc_str} = {term_result}")
        term_counts.append(term_result)
        total_count += term_result

    print("-" * 50)
    
    # Final equation and result
    final_equation = " + ".join(map(str, term_counts))
    print(f"Total number is the sum of these counts:")
    print(f"{final_equation} = {total_count}")

solve()