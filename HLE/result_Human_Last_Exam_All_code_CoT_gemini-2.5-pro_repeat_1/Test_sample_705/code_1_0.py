import itertools

def get_d_operator(states):
    """
    Applies the D operator to a set of states.
    D(S) = union over s in S of {v1, v2, ..., vn}
    """
    decomposed_set = set()
    for state in states:
        decomposed_set.update(state)
    return decomposed_set

def get_c_operator(value_set, all_v):
    """
    Applies the C operator to a set of values.
    C(D) = {(v1, ..., vn) | vk in (D intersect Vk) or (D intersect Vk is empty and vk in Vk)}
    """
    component_options = []
    for k, vk_set in enumerate(all_v):
        intersection = value_set.intersection(vk_set)
        if intersection:
            component_options.append(list(intersection))
        else:
            component_options.append(list(vk_set))

    # itertools.product computes the Cartesian product
    return set(itertools.product(*component_options))

def main():
    """
    Main function to test the claims, specifically to find a counterexample for C.
    """
    # 1. Define the state space components
    V1 = {'v1_a', 'v1_b'}
    V2 = {'v2_x', 'v2_y'}
    ALL_V = [V1, V2]
    S = set(itertools.product(*ALL_V)) # The full state space

    # 2. Define a non-identity function f
    s0 = ('v1_a', 'v2_x')
    s1 = ('v1_b', 'v2_x') # s1 is "close" to s0

    def f(s):
        if s == s0:
            return s1
        return s # Identity for all other states

    # Check if f is identity
    is_identity = all(s == f(s) for s in S)
    print(f"Is the function f an identity function? {is_identity}\n")

    N = 2 # Number of simulation steps

    # 3. Perform ordinary simulation
    ord_sim_states = set()
    current_s = s0
    for _ in range(N + 1):
        ord_sim_states.add(current_s)
        current_s = f(current_s)
    
    print("--- Ordinary Simulation ---")
    print(f"Set of states visited: {sorted(list(ord_sim_states))}\n")

    # 4. Perform relaxed simulation
    sigma = get_d_operator({s0})
    print("--- Relaxed Simulation ---")
    print(f"sigma_0 = {sorted(list(sigma))}")

    for i in range(N):
        # Re-compose states from sigma_i
        recomposed_states = get_c_operator(sigma, ALL_V)
        # Apply f to all re-composed states and decompose the results
        new_values = get_d_operator({f(s) for s in recomposed_states})
        # Update sigma
        sigma.update(new_values)
        print(f"sigma_{i+1} = {sorted(list(sigma))}")

    sigma_N = sigma
    print(f"\nFinal relaxed set sigma_{N}: {sorted(list(sigma_N))}")

    # 5. Apply C to the result of the relaxed simulation
    c_of_sigma_N = get_c_operator(sigma_N, ALL_V)
    print(f"C(sigma_{N}): {sorted(list(c_of_sigma_N))}\n")

    # 6. Compare the results
    print("--- Conclusion for Claim C ---")
    print("Claim C: We can obtain the exactly same result of the ordinary simulation "
          "by applying C to the result of the relaxed simulation if and only if f is identity.")
    
    results_are_equal = (ord_sim_states == c_of_sigma_N)
    print(f"\nDoes ordinary_result == C(relaxed_result)? {results_are_equal}")

    if results_are_equal and not is_identity:
        print("The equality holds, but the function f was NOT identity.")
        print("This provides a counterexample to the 'only if' part of Claim C.")
        print("Therefore, Claim C is FALSE.")
    else:
        print("Could not find a counterexample with this setup.")

if __name__ == "__main__":
    main()