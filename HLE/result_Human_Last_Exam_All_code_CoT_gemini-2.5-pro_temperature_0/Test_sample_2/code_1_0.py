import math
from collections import defaultdict

def parse_group(s):
    """Parses a string like 'Z+Z2' into a dict {'Z': 1, 2: 1}."""
    if s == "0":
        return defaultdict(int)
    parts = s.split('+')
    group = defaultdict(int)
    for part in parts:
        part = part.strip()
        if part == 'Z':
            group['Z'] += 1
        elif part.startswith('Z'):
            n = int(part[1:])
            group[n] += 1
    return group

def format_group(group):
    """Formats a group dict into a string, e.g., {'Z':1, 2:2} -> 'Z + (Z_2)^2'."""
    if not group or all(v == 0 for v in group.values()):
        return "0"
    
    terms = []
    if group.get('Z', 0) > 0:
        count = group['Z']
        if count == 1:
            terms.append("Z")
        else:
            terms.append(f"Z^{count}")
            
    torsion_keys = sorted([k for k in group.keys() if k != 'Z'])
    for n in torsion_keys:
        count = group[n]
        if count > 0:
            if count == 1:
                terms.append(f"Z_{n}")
            else:
                terms.append(f"(Z_{n})^{count}")
            
    return " + ".join(terms)

def add_groups(g1, g2):
    """Computes the direct sum of two abelian groups."""
    res = g1.copy()
    for k, v in g2.items():
        res[k] += v
    return res

def tensor_groups(g1, g2):
    """Computes the tensor product of two abelian groups."""
    res = defaultdict(int)
    for n1, c1 in g1.items():
        for n2, c2 in g2.items():
            count = c1 * c2
            if n1 == 'Z': # Z tensor G = G
                res[n2] += count
            elif n2 == 'Z': # G tensor Z = G
                res[n1] += count
            else: # Z_n tensor Z_m = Z_gcd(n,m)
                g = math.gcd(n1, n2)
                if g > 1:
                    res[g] += count
    return res

def tor_groups(g1, g2):
    """Computes the Tor product of two abelian groups."""
    res = defaultdict(int)
    for n1, c1 in g1.items():
        for n2, c2 in g2.items():
            if n1 == 'Z' or n2 == 'Z': # Tor(Z, G) = 0
                continue
            # Tor(Z_n, Z_m) = Z_gcd(n,m)
            count = c1 * c2
            g = math.gcd(n1, n2)
            if g > 1:
                res[g] += count
    return res

def main():
    """
    Main function to compute the reduced 12-th Spin bordism of BG2.
    """
    # Known Spin bordism groups of a point
    spin_bordism_point = {
        0: "Z", 1: "Z2", 2: "Z2", 3: "0", 4: "Z", 5: "0", 6: "0", 7: "0",
        8: "Z+Z", 9: "Z2+Z2", 10: "Z2", 11: "0", 12: "Z"
    }

    # Known integral homology groups of BG2
    homology_bg2 = {
        0: "Z", 1: "0", 2: "0", 3: "0", 4: "Z", 5: "Z2", 6: "0", 7: "Z2",
        8: "Z+Z2", 9: "Z2+Z2", 10: "Z2", 11: "Z2", 12: "Z+Z2"
    }

    total_group = defaultdict(int)
    equation_terms = []

    # Sum from p=1 to 12 for the reduced group
    for p in range(1, 13):
        q = 12 - p
        
        omega_q_str = spin_bordism_point.get(q, "0")
        if omega_q_str == "0":
            continue
            
        hp_str = homology_bg2.get(p, "0")
        hp_minus_1_str = homology_bg2.get(p - 1, "0")
        
        if hp_str == "0" and hp_minus_1_str == "0":
            continue

        omega_q = parse_group(omega_q_str)
        hp = parse_group(hp_str)
        hp_minus_1 = parse_group(hp_minus_1_str)
        
        # H_p(BG2; Omega_q) = (H_p(BG2) tensor Omega_q) + Tor(H_{p-1}(BG2), Omega_q)
        tensor_part = tensor_groups(hp, omega_q)
        tor_part = tor_groups(hp_minus_1, omega_q)
        
        term_group = add_groups(tensor_part, tor_part)
        
        if any(v > 0 for v in term_group.values()):
            term_str = format_group(term_group)
            equation_terms.append((f"H_{p}(BG_2; \u03A9_{q}^{{spin}})", term_str))
            total_group = add_groups(total_group, term_group)

    # Format the final output
    final_result_str = format_group(total_group)

    print("The reduced 12-th dimensional Spin bordism of BG2 is calculated as the sum:")
    print("\u03A9\u0303_{12}^{spin}(BG_2) = \u2211_{p=1}^{12} H_p(BG_2; \u03A9_{12-p}^{spin})")
    print("\nThe non-zero terms contributing to this sum are:")

    term_expressions = [name for name, result in equation_terms]
    result_expressions = [result.replace("Z", "\u2124") for name, result in equation_terms]
    
    print(" + ".join(term_expressions))
    print("= " + " + ".join(f"({res})" for res in result_expressions))
    
    final_result_formatted = final_result_str.replace("Z", "\u2124")
    print(f"\n= {final_result_formatted}")

if __name__ == "__main__":
    main()