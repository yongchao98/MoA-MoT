import itertools

def get_C_op(D_set, V_map):
    """
    Computes the Composition operator C.
    This is equivalent to the cartesian product of the intersections of D_set with each V_k.
    """
    # Assuming V_k are disjoint, we can find which V_k each value belongs to.
    value_to_k = {v: k for k, V in V_map.items() for v in V}

    # Group the elements of D_set by their V_k affiliation
    # D_k = D_set intersect V_k
    D_k_groups = {k: [] for k in V_map}
    for v in D_set:
        if v in value_to_k:
            D_k_groups[value_to_k[v]].append(v)
    
    # Rule 1: If a group is empty, fill with all possibilities from that V_k
    for k in V_map:
        if not D_k_groups[k]:
            D_k_groups[k] = list(V_map[k])
            
    # Get the component lists in order
    component_lists = [D_k_groups[k] for k in sorted(V_map.keys())]
    
    # The result is the cartesian product
    # Using frozenset for states so they can be added to a set
    C_set = {tuple(p) for p in itertools.product(*component_lists)}
    return C_set


def get_D_op(state):
    """Computes the Decomposition operator D for a single state."""
    return set(state)

# --- Setup a problem instance ---
V1 = {'a1', 'b1'}
V2 = {'a2', 'b2'}
V_map = {1: V1, 2: V2}
S = get_C_op(set(), V_map) # All possible states

# --- Case 1: f is the identity function ---
f_identity = lambda s: s
s0 = ('a1', 'a2')
s1 = f_identity(s0)

print("--- Case 1: f is identity ---")
print(f"s0 = {s0}")
print(f"s1 = f(s0) = {s1}")

# Calculate sigma_1 = D(s0) U D(s1)
sigma_1 = get_D_op(s0).union(get_D_op(s1))
print(f"sigma_1 = D(s0) U D(s1) = {sigma_1}")

# Apply C to sigma_1
C_sigma_1 = get_C_op(sigma_1, V_map)
print(f"C(sigma_1) = {C_sigma_1}")
print(f"Does C(sigma_1) == {{s1}}? {C_sigma_1 == {s1}}")
print("-" * 20)


# --- Case 2: f is NOT the identity function ---
# Example f that swaps the first component
def f_non_identity(s):
    v1, v2 = s
    new_v1 = 'b1' if v1 == 'a1' else 'a1'
    return (new_v1, v2)

s0 = ('a1', 'a2')
s1 = f_non_identity(s0)

print("--- Case 2: f is NOT identity ---")
print(f"s0 = {s0}")
print(f"s1 = f(s0) = {s1}")

# Calculate sigma_1 = D(s0) U D(s1)
sigma_1 = get_D_op(s0).union(get_D_op(s1))
print(f"sigma_1 = D(s0) U D(s1) = {sigma_1}")

# Apply C to sigma_1
C_sigma_1 = get_C_op(sigma_1, V_map)
print(f"C(sigma_1) = {C_sigma_1}")
print(f"Does C(sigma_1) == {{s1}}? {C_sigma_1 == {s1}}")
print("Note: The condition fails as C(sigma_1) contains both s0 and s1.")
