import numpy as np
from itertools import combinations

def get_roots(n, root_type='D'):
    """Generates a set of root vectors for a given type and rank."""
    roots = set()
    if n < 2 and (root_type == 'D' or root_type == 'A'):
        return roots
    if root_type == 'D':
        # Roots are vectors of form +/- e_i +/- e_j, for D_n
        for i in range(n):
            for j in range(i + 1, n):
                for sign1 in [-1, 1]:
                    for sign2 in [-1, 1]:
                        root = np.zeros(n, dtype=int)
                        root[i] = sign1
                        root[j] = sign2
                        roots.add(tuple(root))
    elif root_type == 'A':
        # Roots are vectors of form e_i - e_j, for A_{n-1}
        for i in range(n):
            for j in range(n):
                if i == j: continue
                root = np.zeros(n, dtype=int)
                root[i] = 1
                root[j] = -1
                roots.add(tuple(root))
    return roots

def find_components(roots):
    """Partitions a set of roots into connected components."""
    root_list = [np.array(r) for r in roots]
    num_roots = len(root_list)
    if num_roots == 0:
        return {}

    adj = {i: [] for i in range(num_roots)}
    for i in range(num_roots):
        for j in range(i + 1, num_roots):
            if np.dot(root_list[i], root_list[j]) != 0:
                adj[i].append(j)
                adj[j].append(i)

    visited = [False] * num_roots
    components = {}
    comp_id = 0
    for i in range(num_roots):
        if not visited[i]:
            comp_id += 1
            q = [i]
            visited[i] = True
            component_roots = {tuple(root_list[i])}
            while q:
                u_idx = q.pop(0)
                for v_idx in adj[u_idx]:
                    if not visited[v_idx]:
                        visited[v_idx] = True
                        q.append(v_idx)
                        component_roots.add(tuple(root_list[v_idx]))
            components[comp_id] = component_roots
    return components

def identify_system_type(roots):
    """Identifies the type of a root system based on its structure."""
    num_roots = len(roots)
    if num_roots == 0: return "Empty"
    
    support = set()
    for r in roots:
        for i, coord in enumerate(r):
            if coord != 0: support.add(i)
    
    k = len(support)
    if k==0: return "Empty"
    
    # Remap indices to 0..k-1 to generate canonical roots for comparison
    map_idx = {orig_idx: new_idx for new_idx, orig_idx in enumerate(sorted(list(support)))}
    
    remapped_roots = set()
    for r in roots:
        new_r = np.zeros(k, dtype=int)
        for i, coord in enumerate(r):
            if coord != 0:
                new_r[map_idx[i]] = coord
        remapped_roots.add(tuple(new_r))
        
    # Check for D_k
    if remapped_roots == get_roots(k, 'D'):
        return f"D_{k}"
        
    # Check for A_{k-1}
    if remapped_roots == get_roots(k, 'A'):
        return f"A_{k-1}"
    
    return f"Unknown (rank {k}, {num_roots} roots)"

# --- Problem 1 ---
n1 = 12
u1 = np.ones(n1, dtype=int)
d1 = np.dot(u1, u1)
m_roots1 = {r for r in get_roots(n1, 'D') if np.dot(u1, np.array(r)) % d1 == 0}
sys_type1 = identify_system_type(m_roots1)
ans1 = "Yes" if sys_type1 == "A_11" else "No"

# --- Problem 2 ---
n2 = 15
u2 = np.zeros(n2, dtype=int)
u2[7:] = 1 # 7 zeros, 8 ones
d2 = np.dot(u2, u2)
m_roots2 = {r for r in get_roots(n2, 'D') if np.dot(u2, np.array(r)) % d2 == 0}
components2 = find_components(m_roots2)
has_d7 = False
for cid, c_roots in components2.items():
    if identify_system_type(c_roots) == "D_7":
        has_d7 = True
        break
ans2 = "yes" if has_d7 else "no"

# --- Problem 3 ---
n3 = 18
d3 = 5
possible3 = False
# Primitive vectors u in Z^18 with u.u=5 are permutations of (2,1,0...) or (1,1,1,1,1,0...)
# The structure of the root system only depends on the set of values in u, not their positions.
u_vectors_base = [
    np.array([2,1] + [0]*(n3-2)),
    np.array([1]*5 + [0]*(n3-5))
]
for u3 in u_vectors_base:
    m_roots3 = {r for r in get_roots(n3, 'D') if np.dot(u3, np.array(r)) % d3 == 0}
    components3 = find_components(m_roots3)
    d_count = 0
    for cid, c_roots in components3.items():
        if identify_system_type(c_roots).startswith("D_"):
            d_count += 1
    if d_count > 1:
        possible3 = True
        break
ans3 = "no" if not possible3 else "yes"

# Print the final answer in the requested format
final_answer_str = f"(a) {ans1}; (b) {ans2}; (c) {ans3}"
print(final_answer_str)
print(f"<<<{final_answer_str}>>>")