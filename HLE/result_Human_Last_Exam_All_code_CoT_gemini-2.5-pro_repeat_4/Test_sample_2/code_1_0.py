# A simple representation for finitely generated abelian groups of the form Z^r + (Z_2)^t
# as a tuple (r, t).
AbelianGroup = tuple[int, int]

def add_groups(g1: AbelianGroup, g2: AbelianGroup) -> AbelianGroup:
    """Adds two abelian groups."""
    return (g1[0] + g2[0], g1[1] + g2[1])

def tensor_product(g1: AbelianGroup, g2: AbelianGroup) -> AbelianGroup:
    """Computes the tensor product of two abelian groups."""
    r1, t1 = g1
    r2, t2 = g2
    # (Z^r1 + Z2^t1) @ (Z^r2 + Z2^t2) = 
    # Z^(r1*r2) + Z2^(r1*t2) + Z2^(t1*r2) + Z2^(t1*t2)
    free_rank = r1 * r2
    two_torsion_rank = r1 * t2 + t1 * r2 + t1 * t2
    return (free_rank, two_torsion_rank)

def tor_group(g1: AbelianGroup, g2: AbelianGroup) -> AbelianGroup:
    """Computes Tor_1^Z(g1, g2)."""
    r1, t1 = g1
    r2, t2 = g2
    # Tor(Z, A) = 0, Tor(Z2, Z2) = Z2
    two_torsion_rank = t1 * t2
    return (0, two_torsion_rank)

def hom_to_Z(g: AbelianGroup) -> AbelianGroup:
    """Computes Hom(g, Z)."""
    # Hom(Z, Z) = Z, Hom(Z2, Z) = 0
    return (g[0], 0)

def ext_from_Z(g: AbelianGroup) -> AbelianGroup:
    """Computes Ext^1_Z(g, Z)."""
    # Ext(Z, Z) = 0, Ext(Z2, Z) = Z2
    return (g[1], 0)

def group_to_string(g: AbelianGroup) -> str:
    """Converts a group to a string representation."""
    r, t = g
    if r == 0 and t == 0:
        return "0"
    
    parts = []
    if r > 0:
        parts.append(f"Z^{r}" if r > 1 else "Z")
    if t > 0:
        parts.append(f"(Z_2)^{t}" if t > 1 else "Z_2")
    
    return " + ".join(parts)

def compute_bordism():
    """
    Computes the reduced 12-th Spin bordism group of BG2.
    """
    # Known Spin bordism groups Omega_q^{Spin} = Z^r + (Z_2)^t
    omega_spin = {
        0: (1, 0), 1: (0, 1), 2: (0, 1), 3: (0, 0), 4: (1, 0),
        5: (0, 0), 6: (0, 0), 7: (0, 0), 8: (2, 0), 9: (0, 2),
        10: (0, 3), 11: (0, 0), 12: (3, 0),
    }

    # Known integral cohomology groups H^n(BG2; Z)
    h_cohomology = {
        0: (1, 0), 4: (1, 0), 6: (0, 1), 7: (0, 1), 8: (1, 0),
        9: (0, 2), 10: (0, 2), 11: (0, 2), 12: (2, 0), 13: (0, 2),
    }

    # Step 1: Compute integral homology H_n(BG2; Z) from cohomology using UCT
    # H_n(X;Z) = Hom(H^n(X;Z), Z) + Ext(H^{n+1}(X;Z), Z)
    h_homology = {}
    for n in range(14):
        hn = h_cohomology.get(n, (0, 0))
        hn_plus_1 = h_cohomology.get(n + 1, (0, 0))
        h_homology[n] = add_groups(hom_to_Z(hn), ext_from_Z(hn_plus_1))

    # Step 2: Compute the E^2_{p,q} terms for the reduced group (p > 0)
    # E^2_{p,q} = H_p(BG2; Omega_q) = (H_p(Z) @ Omega_q) + Tor(H_{p-1}(Z), Omega_q)
    total_group = (0, 0)
    equation_terms = []
    
    print("Calculating the filtration quotients E^2_{p, 12-p} for p > 0:")
    for p in range(1, 13):
        q = 12 - p
        omega_q = omega_spin.get(q, (0, 0))
        
        if omega_q == (0, 0):
            continue

        hp = h_homology.get(p, (0, 0))
        hp_minus_1 = h_homology.get(p - 1, (0, 0))
        
        term1 = tensor_product(hp, omega_q)
        term2 = tor_group(hp_minus_1, omega_q)
        e2_pq = add_groups(term1, term2)
        
        if e2_pq != (0, 0):
            print(f"p={p}, q={q}: E^2_{{{p},{q}}} = {group_to_string(e2_pq)}")
            total_group = add_groups(total_group, e2_pq)
            term_str = group_to_string(e2_pq)
            if '+' in term_str:
                term_str = f"({term_str})"
            equation_terms.append(term_str)

    # Step 3: Print the final result
    print("\nThe reduced 12-th Spin bordism group is the direct sum of these terms.")
    final_equation = " + ".join(equation_terms)
    print(f"Final Equation: {final_equation}")
    
    final_result_str = group_to_string(total_group)
    print(f"\nTotal Group: {final_result_str}")
    return final_result_str

if __name__ == '__main__':
    final_answer = compute_bordism()
    # The final answer is requested in a specific format at the end.
    # print(f"\n<<<Z^{total_group[0]} + (Z_2)^{total_group[1]}>>>")
    # This cannot be evaluated in the final block, so I will hardcode the result.

compute_bordism()