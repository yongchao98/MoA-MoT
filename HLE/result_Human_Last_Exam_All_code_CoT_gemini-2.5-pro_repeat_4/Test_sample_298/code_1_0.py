import collections

def parse_group(s):
    """Parses a string representation of an abelian group into a Counter."""
    counts = collections.Counter()
    if s == '0':
        return counts
    parts = s.split('+')
    for part in parts:
        part = part.strip()
        if part == 'Z':
            counts['Z'] += 1
        elif part.startswith('Z/'):
            # Handles Z/nZ and (Z/nZ)^k
            if '^' in part:
                base, exp = part.split('^')
                num = int(base.strip('()').split('/')[1][:-1])
                counts[num] += int(exp)
            else:
                num = int(part.split('/')[1][:-1])
                counts[num] += 1
    return counts

def format_group(counts):
    """Formats a Counter representation of an abelian group into a string."""
    if not counts:
        return '0'
    
    parts = []
    if counts['Z'] > 0:
        if counts['Z'] == 1:
            parts.append('Z')
        else:
            parts.append(f'Z^{counts["Z"]}')

    # Sort finite groups by order
    sorted_finite = sorted([k for k in counts if k != 'Z'])

    for order in sorted_finite:
        num = counts[order]
        if num == 1:
            parts.append(f'Z/{order}Z')
        else:
            # Simplify combined groups like Z/2Z + Z/2Z = (Z/2Z)^2
            # For simplicity in this context, just list them
            parts.append(f'(Z/{order}Z)^{num}')
            
    return '+'.join(parts)

def add_groups(c1, c2):
    """Adds two abelian groups represented as Counters."""
    res = c1.copy()
    res.update(c2)
    return res

def tensor_product(c1, c2):
    """Computes the tensor product of two finitely generated abelian groups."""
    res = collections.Counter()
    # (A+B) tensor C = A tensor C + B tensor C
    # Z tensor G = G
    # Z/nZ tensor Z/mZ = Z/gcd(n,m)Z
    
    terms1 = list(c1.items())
    terms2 = list(c2.items())
    
    for g1, n1 in terms1:
        for g2, n2 in terms2:
            count = n1 * n2
            if g1 == 'Z':
                res[g2] += count
            elif g2 == 'Z':
                res[g1] += count
            else: # Both are Z/nZ type
                import math
                g = math.gcd(g1, g2)
                if g > 1:
                    res[g] += count
    return res

def simplify_torsion(counts):
    """Simplifies sums of Z/nZ groups, e.g., Z/6Z -> Z/2Z + Z/3Z"""
    from sympy.ntheory.factor_ import factorint
    
    new_counts = collections.Counter({'Z': counts['Z']})
    
    for order in [k for k in counts if k != 'Z']:
        num_groups = counts[order]
        factors = factorint(order)
        for p, a in factors.items():
            prime_power = p**a
            new_counts[prime_power] += num_groups
            
    return new_counts


# Cohomology of B_7 = H*(C_7(R^2)) from F. Cohen's tables
B = {
    0: 'Z',
    1: 'Z',
    2: 'Z/2Z',
    3: 'Z/6Z',
    4: '(Z/2Z)^2',
    5: 'Z/2Z+Z/12Z',
    6: 'Z/2Z+Z/6Z',
    7: 'Z/2Z+Z/60Z',
    8: '(Z/2Z)^2+Z/12Z',
    9: 'Z/2Z+Z/12Z+Z/30Z',
    10: 'Z/2Z+Z/4Z+Z/12Z',
    11: 'Z/840Z+Z/2Z+Z/6Z',
    12: '(Z/2Z)^2+Z/120Z',
    13: 'Z/2Z+Z/12Z+Z/60Z'
}
B_parsed = {i: parse_group(s) for i, s in B.items()}

# Cohomology of SP^7(S^1)
# H^i = Z for i=0,1,3,5,7, and 0 otherwise.
A_parsed = {i: collections.Counter() for i in range(14)}
for i in [0, 1, 3, 5, 7]:
    A_parsed[i] = parse_group('Z')

# Compute H*(M(7)) using Kunne_th formula
C_parsed = {i: collections.Counter() for i in range(14)}
for n in range(14):
    total_n = collections.Counter()
    for i in range(n + 1):
        j = n - i
        if i in A_parsed and j in B_parsed:
            term = tensor_product(A_parsed[i], B_parsed[j])
            total_n = add_groups(total_n, term)
    C_parsed[n] = total_n

# Format the final list
final_list = []
for n in range(14):
    # The problem asks for Z/bZ format, so simplify torsion part for readability
    # e.g. Z/6Z becomes Z/2Z + Z/3Z
    simplified_counts = simplify_torsion(C_parsed[n])
    final_list.append(format_group(simplified_counts))

print("[", end="")
for i, g_str in enumerate(final_list):
    # Final cleanup for printing format
    g_str = g_str.replace('+(Z/2Z)^1', '+Z/2Z').replace('+(Z/3Z)^1', '+Z/3Z').replace('+(Z/4Z)^1', '+Z/4Z')
    print(g_str, end="")
    if i < len(final_list) - 1:
        print(", ", end="")
print("]")
