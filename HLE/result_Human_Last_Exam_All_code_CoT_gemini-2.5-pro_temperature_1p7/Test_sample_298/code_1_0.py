from collections import defaultdict

def to_string(group):
    """Converts a group structure to the specified string format."""
    if not group['free'] and not group['torsion']:
        return "0"
    
    parts = []
    if group['free'] > 0:
        if group['free'] == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{group['free']}")
    
    torsion_parts = []
    # Sort torsion components for consistent output
    for n in sorted(group['torsion'].keys()):
        count = group['torsion'][n]
        if count == 1:
            torsion_parts.append(f"Z/{n}Z")
        else:
            torsion_parts.append(f"(Z/{n}Z)^{count}")
            
    parts.extend(torsion_parts)
    return "+".join(parts)

def add_groups(g1, g2):
    """Adds two abelian group structures."""
    res = defaultdict(int)
    res['free'] = g1.get('free', 0) + g2.get('free', 0)
    torsion_res = defaultdict(int)
    for n, c in g1.get('torsion', {}).items():
        torsion_res[n] += c
    for n, c in g2.get('torsion', {}).items():
        torsion_res[n] += c
    res['torsion'] = dict(torsion_res)
    return dict(res)

def ext(hom_group):
    """Computes Ext(H, Z)."""
    res = defaultdict(int)
    res['torsion'] = defaultdict(int)
    torsion_part = hom_group.get('torsion', {})
    for n, count in torsion_part.items():
        res['torsion'][n] += count
    res['torsion'] = dict(res['torsion'])
    return dict(res)

def hom(hom_group):
    """Computes Hom(H, Z)."""
    res = defaultdict(int)
    res['free'] = hom_group.get('free', 0)
    return dict(res)

def parse_group(s):
    """Parses a string representation of a group into a dict structure."""
    if s == "0":
        return {}
    
    g = defaultdict(int)
    g['torsion'] = defaultdict(int)
    
    parts = s.split('+')
    for part in parts:
        if part == "Z":
            g['free'] += 1
        elif part.startswith("Z^"):
            g['free'] += int(part.split('^')[1])
        elif part.startswith("Z/"):
            n = int(part.split('/')[1][:-1])
            g['torsion'][n] += 1
    g['torsion'] = dict(g['torsion'])
    return dict(g)

def main():
    k = 7
    # Homology of Braid groups with sign-twisted coefficients H_p(B_j, Z_sgn)
    # h[j][p] stores H_p(B_j, Z_sgn)
    h_str = {
        1: {0: "Z"},
        2: {0: "Z/2Z"},
        3: {1: "Z/2Z", 3: "Z/6Z"},
        4: {0: "Z/2Z", 2: "Z/2Z", 4: "Z/10Z"},
        5: {1: "Z/2Z", 3: "Z/2Z+Z/2Z", 5: "Z/6Z+Z/6Z"},
        6: {0: "Z/2Z+Z/2Z", 2: "Z/2Z+Z/2Z+Z/2Z", 4: "Z/2Z+Z/2Z", 6: "Z/14Z+Z/30Z"},
        7: {1: "Z/3Z+Z/2Z", 3: "Z/2Z", 5: "Z/6Z+Z/10Z"}
    }
    
    h = defaultdict(lambda: defaultdict(dict))
    for j, data in h_str.items():
        for p, s in data.items():
            h[j][p] = parse_group(s)

    max_dim = 0
    for j in range(1, k + 1):
        for p in h[j]:
            max_dim = max(max_dim, p + j)

    homology = defaultdict(dict)
    for i in range(max_dim + 1):
        current_homology = {}
        for j in range(1, k + 1):
            p = i - j
            if p >= 0 and p in h[j]:
                current_homology = add_groups(current_homology, h[j][p])
        homology[i] = current_homology
    
    cohomology = []
    for i in range(max_dim + 2):
        h_i = homology.get(i, {})
        h_i_minus_1 = homology.get(i - 1, {})
        
        term1 = hom(h_i)
        term2 = ext(h_i_minus_1)
        
        coh_i = add_groups(term1, term2)
        cohomology.append(coh_i)

    # Trim trailing zero groups
    last_nonzero = -1
    for i in range(len(cohomology) - 1, -1, -1):
        if to_string(cohomology[i]) != "0":
            last_nonzero = i
            break
            
    final_cohomology_list = [to_string(g) for g in cohomology[:last_nonzero + 1]]
    
    print(final_cohomology_list)

main()