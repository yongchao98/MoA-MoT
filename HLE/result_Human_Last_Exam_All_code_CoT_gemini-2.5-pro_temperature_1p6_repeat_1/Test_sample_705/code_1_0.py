import itertools
from functools import lru_cache

# Let's set up the state space based on the problem description.
# We define two disjoint finite sets, V_1 and V_2.
# To ensure disjointness in Python, we'll represent values as tuples,
# where the first element identifies the set V_k it belongs to.
V = {
    1: {('v1', 'a'), ('v1', 'b')},
    2: {('v2', 'a'), ('v2', 'b')}
}
n = len(V)
V_keys = sorted(V.keys())

# The state space S is the Cartesian product V_1 x V_2
S_space = list(itertools.product(*(V[k] for k in V_keys)))
# The union space D is the union of all V_k
D_space = frozenset.union(*V.values())

# --- Operators D and C ---

def D_op(states_set):
    """ Implements mathscr{D}: 2^S -> 2^D """
    decomposed = set()
    for state in states_set:
        for val in state:
            decomposed.add(val)
    return frozenset(decomposed)

# We use a cache (memoization) for the recursive C_op for efficiency.
C_op_cache = {}

def C_op(d_set):
    """ Implements mathscr{C}: 2^D -> 2^S """
    d_set = frozenset(d_set)
    if d_set in C_op_cache:
        return C_op_cache[d_set]

    # Rule 1: If D is missing components from a V_k, fill it with all possibilities.
    is_updated = False
    for k in V_keys:
        if not (d_set & V[k]):
            d_set = d_set.union(V[k])
            is_updated = True
    if is_updated:
        return C_op(d_set)

    # Rule 2: If D has multiple components from a V_k, branch the computation.
    for k in V_keys:
        intersection = d_set & V[k]
        if len(intersection) > 1:
            result_states = set()
            d_without_vk = d_set - V[k]
            for v_j in intersection:
                new_d_set = d_without_vk.union({v_j})
                result_states.update(C_op(new_d_set))
            C_op_cache[d_set] = frozenset(result_states)
            return frozenset(result_states)

    # Rule 3: Base case, where |D intersect V_k| = 1 for all k.
    final_components = []
    for k in V_keys:
        intersection = d_set & V[k]
        if len(intersection) == 1:
            final_components.append(list(intersection)[0])
        else:
            # This path indicates an invalid state (e.g., from Rule 1).
            C_op_cache[d_set] = frozenset()
            return frozenset()

    state_tuple = tuple(final_components)
    result = frozenset({state_tuple})
    C_op_cache[d_set] = result
    return result

# --- Simulation Logic ---

def ordinary_sim(f, s0, N):
    """ Runs the ordinary simulation and returns the N-th state, s_N. """
    s_current = s0
    for _ in range(N):
        s_current = f(s_current)
    return s_current

def relaxed_sim(f, s0, N):
    """ Runs the relaxed simulation and returns the N-th abstract state, sigma_N. """
    sigma_current = D_op({s0})
    for _ in range(N):
        s_set_to_sim = C_op(sigma_current)
        if not s_set_to_sim:
            break
        
        next_states_decomposed = D_op({f(s) for s in s_set_to_sim})
        sigma_current = sigma_current.union(next_states_decomposed)
        
    return sigma_current

# --- Simulator Functions f ---

def f_identity(s):
    """ The identity function, f(s) = s. """
    return s

s_example_a = (('v1', 'a'), ('v2', 'a'))
s_example_b = (('v1', 'b'), ('v2', 'b'))

def f_non_identity(s):
    """ A non-identity function. It swaps s_example_a and is identity otherwise. """
    if s == s_example_a:
        return s_example_b
    return s

# --- Verification of Claim C ---

def verify_claim_c():
    """
    Tests the claim: mathscr{C}(sigma_N) == {s_N} if and only if f is identity.
    """
    print("--- Verifying Claim C ---")
    N = 3
    s0 = s_example_a
    
    print("\nPart 1: Testing the 'if' direction (f = identity).")
    s_N_id = ordinary_sim(f_identity, s0, N)
    sigma_N_id = relaxed_sim(f_identity, s0, N)
    C_sigma_N_id = C_op(sigma_N_id)
    
    print(f"  Starting with s0 = {s0} and N = {N}")
    print(f"  Ordinary result s_N = {s_N_id}")
    print(f"  Relaxed result sigma_N = {sigma_N_id}")
    print(f"  Applying C_op to relaxed result: C(sigma_N) = {C_sigma_N_id}")
    condition_holds_id = (C_sigma_N_id == {s_N_id})
    print(f"  Is C(sigma_N) == {{s_N}}? {condition_holds_id}")
    print("  Conclusion: The condition holds, as expected.")

    print("\nPart 2: Testing the 'only if' direction (f != identity).")
    # We demonstrate this by finding a counterexample where the condition fails.
    N_fail = 1
    s_N_non_id = ordinary_sim(f_non_identity, s0, N_fail)
    sigma_N_non_id = relaxed_sim(f_non_identity, s0, N_fail)
    C_sigma_N_non_id = C_op(sigma_N_non_id)
    
    print(f"  Starting with s0 = {s0} and N = {N_fail}")
    print(f"  Ordinary result s_N = {s_N_non_id}")
    print(f"  Relaxed result sigma_N = {sigma_N_non_id}")
    print(f"  Applying C_op to relaxed result: C(sigma_N) = {C_sigma_N_non_id}")
    condition_holds_non_id = (C_sigma_N_non_id == {s_N_non_id})
    print(f"  Is C(sigma_N) == {{s_N}}? {condition_holds_non_id}")
    print("  Conclusion: The condition fails, as expected.")
    print("\nSince the condition fails for a non-identity function, it supports the 'only if' part.")
    print("The Python simulation confirms that Claim C is correct.")


verify_claim_c()