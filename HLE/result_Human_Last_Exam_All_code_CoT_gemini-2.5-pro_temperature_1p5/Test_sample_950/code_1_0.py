def solve_grassmannian_torsion():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the Grassmannian Gr(3,5).

    The rank is interpreted as the number of cyclic summands in the torsion subgroup.
    """
    # The integral homology groups H_i(Gr(3,5); Z) for i = 0 to 6.
    # Each group is represented as a tuple (rank_of_free_part, list_of_torsion_coeffs).
    # Based on Korbaš and Lörinc, "On the cohomology of real Grassmann manifolds".
    # H_i = Z^{rank} + Z_{c_1} + Z_{c_2} + ...
    # H_groups[i] represents H_i
    H_groups = [
        (1, []),          # H_0 = Z
        (1, []),          # H_1 = Z
        (2, []),          # H_2 = Z^2
        (2, [2]),         # H_3 = Z^2 + Z_2
        (2, [2]),         # H_4 = Z^2 + Z_2
        (1, [2]),         # H_5 = Z   + Z_2
        (1, []),          # H_6 = Z
        (0, [])           # H_7 = 0
    ]
    
    total_torsion_rank = 0
    
    print("Calculating the torsion part of each cohomology group H^k(Gr(3,5); Z):")
    # H^k has torsion part isomorphic to the torsion part of H_{k-1}.
    # We sum the number of torsion components for k from 1 up to dim+1.
    # The dimension of Gr(3,5) is 3 * (5-3) = 6. So we check up to k=7.
    for k in range(1, len(H_groups)):
        # H_{k-1} is given by H_groups[k-1]
        H_k_minus_1 = H_groups[k-1]
        
        # Tors(H^k) is isomorphic to Tors(H_{k-1})
        torsion_coeffs = H_k_minus_1[1]
        
        # The number of cyclic summands in the torsion group
        torsion_rank_k = len(torsion_coeffs)
        
        if torsion_rank_k > 0:
            torsion_group_str = " + ".join([f"Z_{c}" for c in torsion_coeffs])
            print(f"Tors(H^{k}) is isomorphic to Tors(H_{k-1}) = {torsion_group_str}. Rank = {torsion_rank_k}")
            total_torsion_rank += torsion_rank_k
        else:
             print(f"Tors(H^{k}) is isomorphic to Tors(H_{k-1}) = 0. Rank = 0")


    print("\nThe total rank of the torsion subgroup of the integral cohomology ring is the sum of these ranks.")
    print(f"The calculation is: 1 (from H^4) + 1 (from H^5) + 1 (from H^6)")
    print(f"Total rank = {total_torsion_rank}")

solve_grassmannian_torsion()