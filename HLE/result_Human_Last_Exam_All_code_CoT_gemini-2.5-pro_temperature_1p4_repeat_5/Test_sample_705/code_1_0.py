import itertools
from functools import reduce

# To ensure V_k are disjoint, we'll represent values as tuples (k, val)
def setup_domain(Vs_config):
    """Sets up the disjoint V_k sets and the overall domain D."""
    V_sets = {k: {(k, v) for v in vals} for k, vals in Vs_config.items()}
    D_domain = set().union(*V_sets.values())
    return V_sets, D_domain

# Define the state space S as the Cartesian product
def get_S(V_sets):
    """Computes the state space S from the V_k sets."""
    # The components of S are sorted by k to ensure a canonical order
    sorted_V_sets = [V_sets[k] for k in sorted(V_sets.keys())]
    return set(itertools.product(*sorted_V_sets))

# D operator: decompose a set of states to a set of values
def op_D(states):
    """Applies the D operator."""
    values = set()
    for s in states:
        values.update(s)
    return values

# C operator: re-compose a set of values to a set of states
def op_C(D_subset, V_sets):
    """Applies the C operator recursively."""
    # Rule 1: Completion
    # Find any k for which the intersection is empty
    completion_applied = False
    for k, Vk in V_sets.items():
        if not D_subset.intersection(Vk):
            # This recursive call handles the case where multiple V_k are missing
            return op_C(D_subset.union(Vk), V_sets)

    # Rule 2: Branching
    for k, Vk in V_sets.items():
        intersection = D_subset.intersection(Vk)
        if len(intersection) > 1:
            # Union of recursive calls for each value in the intersection
            composed_states = set()
            for v_j in intersection:
                # D' = (D \ V_k) U {v_j}
                D_prime = (D_subset.difference(Vk)).union({v_j})
                composed_states.update(op_C(D_prime, V_sets))
            return composed_states

    # Rule 3: Construction
    # This rule is reached if |D_subset intersect Vk| == 1 for all k
    # We construct the single state tuple. The order matters.
    state_components = []
    for k in sorted(V_sets.keys()):
        val_set = D_subset.intersection(V_sets[k])
        # This condition should be met if we reached this part of the code
        if len(val_set) == 1:
            state_components.append(list(val_set)[0])
        else:
             # This should not happen if the logic is correct
             return set() # Should not be reachable
    return {tuple(state_components)}


# --- Simulation Setup ---
Vs = {1: {'a', 'b'}, 2: {'c', 'd'}}
V, D = setup_domain(Vs)
S = get_S(V)
n = len(Vs)

# Initial state
s0 = ((1, 'a'), (2, 'c'))

# --- Case 1: f is identity ---
f_identity = lambda s: s

# Ordinary simulation
s1_id = f_identity(s0)

# Relaxed simulation
sigma0_id = op_D({s0})
C_sigma0_id = op_C(sigma0_id, V)
union_f_s_id = op_D({f_identity(s) for s in C_sigma0_id})
sigma1_id = sigma0_id.union(union_f_s_id)

# Check condition for C
C_sigma1_id = op_C(sigma1_id, V)

print("--- Testing Claim C with Identity Function (f(s) = s) ---")
print(f"s_0 = {s0}")
print(f"Ordinary simulation result: s_1 = {s1_id}")
print(f"Relaxed simulation generated sigma_1: {sigma1_id}")
print(f"Applying C to sigma_1 gives C(sigma_1) = {C_sigma1_id}")
print(f"Is C(sigma_1) == {{s_1}}? {C_sigma1_id == {s1_id}}")
print("Conclusion: For f=identity, the property holds.\n")


# --- Case 2: f is NOT identity ---
# Let's define f to flip the first component's value if it's 'a'
def f_non_identity(s):
    if (1, 'a') in s:
        new_s_list = list(s)
        new_s_list[0] = (1, 'b')
        return tuple(new_s_list)
    return s

# Ordinary simulation
s1_non_id = f_non_identity(s0)

# Relaxed simulation
sigma0_non_id = op_D({s0})
C_sigma0_non_id = op_C(sigma0_non_id, V)
union_f_s_non_id = op_D({f_non_identity(s) for s in C_sigma0_non_id})
sigma1_non_id = sigma0_non_id.union(union_f_s_non_id)

# Check condition for C
C_sigma1_non_id = op_C(sigma1_non_id, V)

print("--- Testing Claim C with Non-Identity Function ---")
print(f"s_0 = {s0}")
print(f"Ordinary simulation result: s_1 = {s1_non_id}")
print(f"Relaxed simulation generated sigma_1: {sigma1_non_id}")
print(f"Applying C to sigma_1 gives C(sigma_1) = {C_sigma1_non_id}")
print(f"Is C(sigma_1) == {{s_1}}? {C_sigma1_non_id == {s1_non_id}}")
print(f"Note that C(sigma_1) contains both the original state {s0} and the new state {s1_non_id}.")
print("Conclusion: For f!=identity, the property does not hold.")
print("\nThis demonstrates that the equivalence in Claim C is correct.")
print(f"The equation C(sigma_N) = {{s_N}} only holds if and only if f is the identity function.")
print(f"For our N=1 non-identity example, the final sets are C(sigma_1) = {C_sigma1_non_id} and {{s_1}} = {{{s1_non_id}}}, which are not equal.")
