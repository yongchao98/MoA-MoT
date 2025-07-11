import math

def group_name(free_rank, z2_rank):
    """Generates a string representation of an abelian group."""
    parts = []
    if free_rank > 0:
        if free_rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{free_rank}")
    if z2_rank > 0:
        if z2_rank == 1:
            parts.append("Z_2")
        else:
            parts.append(f"Z_2^{z2_rank}")
    if not parts:
        return "0"
    return " + ".join(parts)

def tensor_product(g1, g2):
    """Computes the tensor product of two abelian groups of the form (free_rank, z2_rank)."""
    g1_free, g1_z2 = g1
    g2_free, g2_z2 = g2
    
    # Z x Z = Z
    free_part = g1_free * g2_free
    
    # Z x Z_2 = Z_2
    # Z_2 x Z = Z_2
    # Z_2 x Z_2 = Z_2
    z2_part = g1_free * g2_z2 + g1_z2 * g2_free + g1_z2 * g2_z2
    
    return (free_part, z2_part)

def main():
    """
    Computes the reduced 12-th dimensional Spin bordism of the classifying space of G2.
    """
    # Integral homology groups H_p(BG2; Z) for p <= 12.
    # Represented as a tuple (rank of free part, rank of Z_2 part).
    # Other torsion is not present in these dimensions.
    H = {
        0: (1, 0), 1: (0, 0), 2: (0, 0), 3: (0, 0), 4: (1, 0),
        5: (0, 1), 6: (0, 0), 7: (0, 1), 8: (1, 0), 9: (0, 2),
        10: (0, 1), 11: (0, 0), 12: (2, 0)
    }

    # Spin bordism groups Omega_q^Spin for q <= 12.
    Omega = {
        0: (1, 0), 1: (0, 1), 2: (0, 1), 3: (0, 0), 4: (1, 0),
        5: (0, 0), 6: (0, 0), 7: (0, 0), 8: (1, 0), 9: (0, 2),
        10: (0, 1), 11: (0, 0), 12: (1, 0)
    }

    print("Computing the reduced 12-th Spin bordism group of BG2: Sum of H_p(BG2) x Omega_{12-p}(pt) for p > 0.")
    print("-" * 80)
    
    total_free_rank = 0
    total_z2_rank = 0
    
    # Store the non-zero terms for the final equation
    equation_terms = []

    # p ranges from 1 to 12 for the reduced group
    for p in range(1, 13):
        q = 12 - p
        
        hp = H.get(p, (0, 0))
        omega_q = Omega.get(q, (0, 0))
        
        # If either group is trivial, the tensor product is trivial.
        if (hp[0] == 0 and hp[1] == 0) or (omega_q[0] == 0 and omega_q[1] == 0):
            continue

        term_group = tensor_product(hp, omega_q)
        
        total_free_rank += term_group[0]
        total_z2_rank += term_group[1]

        hp_name = group_name(hp[0], hp[1])
        omega_q_name = group_name(omega_q[0], omega_q[1])
        term_name = group_name(term_group[0], term_group[1])
        
        print(f"Term p={p}, q={q}: H_{p} x Omega_{q} = ({hp_name}) x ({omega_q_name}) = {term_name}")
        equation_terms.append(term_name)

    print("-" * 80)
    final_group_name = group_name(total_free_rank, total_z2_rank)
    final_equation = " + ".join(equation_terms)
    print(f"The total reduced bordism group is the direct sum of these terms:")
    print(f"Final Sum = {final_equation} = {final_group_name}")

if __name__ == "__main__":
    main()
