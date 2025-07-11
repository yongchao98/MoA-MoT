def add_groups(g1, g2):
    """Adds two abelian groups represented as dictionaries."""
    res = g1.copy()
    for order, count in g2.items():
        res[order] = res.get(order, 0) + count
    return res

def format_group(g):
    """Formats an abelian group dictionary into the specified string format."""
    if not g or all(v == 0 for v in g.values()):
        return "0"
    
    parts = []
    # Handle the free part (Z)
    if g.get(0, 0) > 0:
        parts.extend(["Z"] * g.get(0, 0))
    
    # Handle the torsion parts (Z/bZ)
    # Sort by order for consistent output
    sorted_torsion_keys = sorted([k for k in g.keys() if k != 0])
    for order in sorted_torsion_keys:
        count = g[order]
        parts.extend([f"Z/{order}Z"] * count)
        
    return "+".join(parts)

def main():
    """
    Computes and prints the list of cohomology groups for M(7).
    """
    # Integral cohomology of the braid group B_7, H^i(B_7; Z).
    # We represent an abelian group like Z^a + (Z/n1Z)^b1 + ...
    # as a dictionary {0: a, n1: b1, ...}.
    H_B7 = {
        0: {0: 1},          # H^0(B_7) = Z
        1: {0: 1},          # H^1(B_7) = Z
        2: {2: 1},          # H^2(B_7) = Z/2Z
        3: {2: 1},          # H^3(B_7) = Z/2Z
        4: {2: 2},          # H^4(B_7) = Z/2Z + Z/2Z
        5: {2: 1},          # H^5(B_7) = Z/2Z
        6: {12: 1},         # H^6(B_7) = Z/12Z
    }
    # H^i(B_7) = 0 for i > 6

    H_M7_list = []
    # Cohomology of M(7) is non-zero up to degree 7
    top_degree = 7
    for n in range(top_degree + 1):
        # H^n(M(7)) = H^n(B_7) + H^{n-1}(B_7)
        H_n_B7 = H_B7.get(n, {})
        H_n_minus_1_B7 = H_B7.get(n - 1, {})
        
        H_n_M7 = add_groups(H_n_B7, H_n_minus_1_B7)
        H_M7_list.append(format_group(H_n_M7))
        
    # Format the final output string as requested
    final_output = f"[{', '.join(H_M7_list)}]"
    print(final_output)

if __name__ == "__main__":
    main()
