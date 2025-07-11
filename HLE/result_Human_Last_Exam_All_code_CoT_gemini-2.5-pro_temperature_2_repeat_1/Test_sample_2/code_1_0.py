import math

def gcd(a, b):
    return math.gcd(a, b)

# Represents a finitely generated abelian group as a dictionary
# {0: rank_of_free_part, 2: num_Z2, 3: num_Z3, ...}
def group_add(g1, g2):
    """Adds two abelian groups."""
    res = g1.copy()
    for k, v in g2.items():
        res[k] = res.get(k, 0) + v
    return res

def group_to_string(g):
    """Formats a group dictionary into a readable string."""
    parts = []
    if g.get(0, 0) > 0:
        if g[0] == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{g[0]}")
    
    sorted_torsion_keys = sorted([k for k in g if k > 1])
    for k in sorted_torsion_keys:
        v = g[k]
        if v > 0:
            if v == 1:
                parts.append(f"Z_{k}")
            else:
                parts.append(f"(Z_{k})^{v}")
    
    if not parts:
        return "0"
    return " + ".join(parts)

# Implements tensor product G H where H is Z, Z_n, or a direct sum of these
def tensor_product(g, h):
    """Computes tensor product of two abelian groups."""
    res = {}
    
    # G tensor Z = G
    res = group_add(res, {k: v * h.get(0,0) for k, v in g.items()})

    for h_n, h_v in h.items():
        if h_n == 0: continue # Z part already done
        for _ in range(h_v):
            # G tensor Z_n = G / nG
            temp_res = {}
            # Z tensor Z_n = Z_n
            if g.get(0,0) > 0:
                temp_res[h_n] = temp_res.get(h_n, 0) + g[0]

            # Z_m tensor Z_n = Z_gcd(m,n)
            for g_m, g_v in g.items():
                if g_m == 0: continue
                d = gcd(g_m, h_n)
                if d > 1:
                    temp_res[d] = temp_res.get(d, 0) + g_v
            res = group_add(res, temp_res)

    return res


def tor_fun(g, h):
    """Computes Tor_1(G, H) for our specific groups."""
    res = {}
    
    # Tor(Z, G) = 0, so free part of g contributes nothing
    
    # Tor(Z_n, Z_m) = Z_gcd(n,m)
    for g_n, g_v in g.items():
        if g_n == 0: continue # Free part done
        for h_m, h_v in h.items():
            if h_m == 0: continue # Free part Tor is 0
            
            d = gcd(g_n, h_m)
            if d > 1:
                res[d] = res.get(d, 0) + g_v * h_v
    return res

# Spin Bordism Groups of a Point
# Omega_q^{Spin}
Omega_spin = [
    {0: 1},        # q=0: Z
    {2: 1},        # q=1: Z_2
    {2: 1},        # q=2: Z_2
    {},            # q=3: 0
    {0: 1},        # q=4: Z
    {},            # q=5: 0
    {},            # q=6: 0
    {},            # q=7: 0
    {0: 2},        # q=8: Z + Z
    {2: 2},        # q=9: Z_2 + Z_2
    {2: 1},        # q=10: Z_2
    {},            # q=11: 0
    {0: 1},        # q=12: Z
]

# Reduced Integral Homology of BG2
# H_p(BG2)
H_BG2 = [
    {},            # p=0 (reduced homology)
    {},            # p=1: 0
    {},            # p=2: 0
    {},            # p=3: 0
    {0: 1},        # p=4: Z
    {},            # p=5: 0
    {2: 1},        # p=6: Z_2
    {},            # p=7: 0
    {0: 1, 3: 1},  # p=8: Z + Z_3
    {2: 1},        # p=9: Z_2
    {2: 1},        # p=10: Z_2
    {},            # p=11: 0
    {0: 1, 2: 1},  # p=12: Z + Z_2
]

def compute_spin_bordism_of_bg2():
    """
    Computes the 12th reduced Spin bordism of BG2 using a collapsing AHSS.
    """
    total_degree = 12
    final_group = {}
    
    equation_parts = []
    
    print("Computing the E^2_{p,q} terms of the Atiyah-Hirzebruch Spectral Sequence for p+q=12.")
    print("Reduced Spin bordism grouptilde(Omega)_12^Spin(BG_2) is the direct sum of these terms.\n")
    
    for p in range(1, total_degree + 1):
        q = total_degree - p
        if q < 0 or q >= len(Omega_spin): continue

        Omega_q = Omega_spin[q]
        H_p = H_BG2[p]
        H_p_minus_1 = H_BG2[p-1] if p > 0 else {}
        
        # H_p(X; G) = (H_p(X) tensor G) + Tor(H_{p-1}(X), G)
        tensor_term = tensor_product(H_p, Omega_q)
        tor_term = tor_fun(H_p_minus_1, Omega_q)
        
        e2_pq = group_add(tensor_term, tor_term)
        
        if e2_pq:
            term_str = group_to_string(e2_pq)
            print(f"E^2_{{{p},{q}}} = H_{p}(BG_2; Omega_{q}^{{Spin}}) = {term_str}")
            final_group = group_add(final_group, e2_pq)
            equation_parts.append(term_str)
        
    print("\n------------------------------------------------------------")
    print("The total reduced Spin bordism group is the sum:")
    print(f"tilde(Omega)_{{12}}^{{Spin}}(BG_2) = {' + '.join(equation_parts)}")
    
    final_string = group_to_string(final_group)
    print("\nFinal Result:")
    print(f"tilde(Omega)_{{12}}^{{Spin}}(BG_2) = {final_string}")
    
    # For the final answer extraction
    final_answer_formatted = final_string.replace(" + ", " ")
    print(f"\n<<<{final_answer_formatted}>>>")


compute_spin_bordism_of_bg2()
