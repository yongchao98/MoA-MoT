import collections

def get_spin_bordism_groups():
    """Returns the first few Spin bordism groups of a point."""
    # Omega^{Spin}_q(pt)
    # Using Z for ZZ, Z2 for ZZ/2ZZ
    # List representing [Z_rank, Z2_count, Z3_count, ...]
    omega = {
        0: [1, 0],
        1: [0, 1],
        2: [0, 1],
        3: [0, 0],
        4: [1, 0],
        5: [0, 0],
        6: [0, 0],
        7: [0, 0],
        8: [2, 0],
        9: [0, 2],
        10: [0, 1],
        11: [0, 0],
        12: [1, 0],
    }
    return omega

def get_bsu3_homology_ranks():
    """Returns the ranks of the integral homology groups of BSU(3)."""
    # H_p(BSU(3); Z) is free abelian.
    # H_*(BSU(3); Z) is dual to Z[c2, c3]
    h_ranks = collections.defaultdict(int)
    h_ranks[0] = 1
    # Combinations of c2 (deg 4) and c3 (deg 6)
    # Dim 4: c2
    h_ranks[4] = 1
    # Dim 6: c3
    h_ranks[6] = 1
    # Dim 8: c2^2
    h_ranks[8] = 1
    # Dim 10: c2*c3
    h_ranks[10] = 1
    # Dim 12: c2^3, c3^2
    h_ranks[12] = 2
    return h_ranks

def tensor_product(g1, g2):
    """Computes a simplified tensor product for Z, Z_n."""
    # g = [Z_rank, Z2_count, Z3_count, ...]
    if g1 == [0, 0] or g2 == [0, 0]:
        return [0, 0]
    # Z tensor G = G
    if g1 == [1, 0]: return g2
    if g2 == [1, 0]: return g1
    # Z_n tensor Z_m = Z_gcd(n,m)
    z_rank = 0 # No Z component from tensor of torsion
    z2_count = g1[1] * g2[1]
    # For this problem we only need Z and Z2
    return [z_rank, z2_count]

def compute_E2_term(p, h_ranks, omega):
    """Computes E^2_{p,q} = H_p(BSU(3); Omega_q)."""
    q = 12 - p
    if p not in h_ranks or q not in omega:
        return [0, 0]
    
    h_p = [h_ranks[p], 0]  # H_p(BSU(3)) is free
    omega_q = omega[q]

    # H_p(X; G) = H_p(X) tensor G (since H_{p-1} is free, Tor term vanishes)
    # The groups are Z, Z2
    
    # Z^k tensor [rank, torsion_2] = [k*rank, k*torsion_2]
    z_rank = h_p[0] * omega_q[0]
    z2_count = h_p[0] * omega_q[1]
    
    return [z_rank, z2_count]

def format_group(group_parts):
    """Formats the group structure into a string."""
    parts = []
    if group_parts['Z'] > 0:
        parts.extend(['Z'] * group_parts['Z'])
    if group_parts['Z2'] > 0:
        parts.extend(['Z_2'] * group_parts['Z2'])
    if not parts:
        return "0"
    return " (+) ".join(parts)

def main():
    """Main computation logic."""
    print("Computing the reduced 12-th dimensional Spin bordism of BG_2.")
    print("This is isomorphic to the reduced 12-th Spin bordism of BSU(3).\n")
    
    omega = get_spin_bordism_groups()
    h_ranks = get_bsu3_homology_ranks()
    target_dim = 12

    total_E2_group = collections.defaultdict(int)
    
    print("Calculating the E^2 page of the Atiyah-Hirzebruch spectral sequence for p+q=12, p>0:")
    print("-" * 70)
    print(f"{'p':<5} {'q':<5} {'H_p(BSU(3))':<15} {'Omega_q':<15} {'E^2_{p,q} = H_p(Omega_q)':<20}")
    print("-" * 70)
    
    # Iterate p from 1 to 12 for the reduced group
    for p in range(1, target_dim + 1):
        q = target_dim - p
        e2_term = compute_E2_term(p, h_ranks, omega)
        
        # Formatting for display
        h_p_str = f"Z^{h_ranks[p]}" if h_ranks[p] > 1 else ("Z" if h_ranks[p] == 1 else "0")
        omega_q_str = format_group({'Z': omega.get(q, [0,0])[0], 'Z2': omega.get(q, [0,0])[1]})
        e2_term_str = format_group({'Z': e2_term[0], 'Z2': e2_term[1]})
        
        if e2_term[0] > 0 or e2_term[1] > 0:
            print(f"{p:<5} {q:<5} {h_p_str:<15} {omega_q_str:<15} {e2_term_str:<20}")

        total_E2_group['Z'] += e2_term[0]
        total_E2_group['Z2'] += e2_term[1]
    
    print("-" * 70)
    e2_sum_str = format_group(total_E2_group)
    print(f"The sum of these E^2 terms is: {e2_sum_str}\n")
    
    print("Considering the differentials in the spectral sequence:")
    print("The term E^2_{10,2} = Z_2 is killed by the d_2 differential from E^2_{12,1}.")
    print("Specifically, d_2: H_{12}(BSU(3); Omega_1) -> H_{10}(BSU(3); Omega_2) is surjective.")
    print("Other differentials on this diagonal are zero.")
    
    final_group = total_E2_group.copy()
    final_group['Z2'] -= 1  # The Z_2 term is eliminated
    
    final_group_str = format_group(final_group)
    print("\nAfter accounting for differentials, the resulting group is a free abelian group.")
    print("The final result for the reduced 12-th Spin bordism of BG_2 is:")
    
    # Create the equation Z + Z + ... from the rank
    z_components = ['Z'] * final_group['Z']
    final_equation = " (+) ".join(z_components)
    print(f"Final Answer = {final_equation}")

if __name__ == "__main__":
    main()
