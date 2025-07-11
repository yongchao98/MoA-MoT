import itertools

def demonstrate_claim_C():
    """
    This script demonstrates the logic behind Claim C by testing the
    condition {s_1} = C(sigma_1) for two cases:
    1. f is the identity function (s_1 = s_0)
    2. f is not the identity function (s_1 != s_0)
    """

    # Let's define the setting for n=3
    # V_k = {(k, 0), (k, 1)} are disjoint finite sets
    n = 3
    V = {k: {(k, 0), (k, 1)} for k in range(1, n + 1)}

    # D operator: Decomposes states into a set of their component values
    def D_op(set_of_states):
        result = set()
        for s in set_of_states:
            result.update(s)
        return result

    # C operator: Composes a set of values into a set of states
    def C_op(set_of_values):
        # This implementation handles the composition based on the rules.
        # It finds the values in the input set for each component V_k
        # and then generates the Cartesian product of these value sets.
        component_options = []
        for k in range(1, n + 1):
            options = set_of_values.intersection(V[k])
            # Rule 1 (not needed for this demo): if not options: options = V[k]
            component_options.append(list(options))
        
        # itertools.product handles Rule 2 (branching) and Rule 3 (base case)
        return set(itertools.product(*component_options))

    # --- Test Case 1: f is identity ---
    print("--- Case 1: f is the identity function (s_1 = s_0) ---")
    s_0_case1 = ((1, 0), (2, 0), (3, 0))
    # If f is identity, s_1 = f(s_0) = s_0
    s_1_case1 = s_0_case1

    # sigma_1 = D(s_0) U D(s_1)
    sigma_1_case1 = D_op({s_0_case1, s_1_case1})
    
    # Check the condition: {s_1} = C(sigma_1)
    C_result_case1 = C_op(sigma_1_case1)
    
    # Print the "equation" parts
    print(f"The state s_1 is: {s_1_case1}")
    print(f"The set of values sigma_1 is: {sigma_1_case1}")
    print(f"The composed set C(sigma_1) is: {C_result_case1}")
    print(f"The condition {{s_1}} == C(sigma_1) holds: {C_result_case1 == {s_1_case1}}\n")


    # --- Test Case 2: f is NOT identity ---
    print("--- Case 2: f is NOT the identity function (s_1 != s_0) ---")
    s_0_case2 = ((1, 0), (2, 0), (3, 0))
    # Let f map s_0 to a different state s_1
    s_1_case2 = ((1, 1), (2, 0), (3, 0)) # Differs in the first component

    # sigma_1 = D(s_0) U D(s_1)
    sigma_1_case2 = D_op({s_0_case2, s_1_case2})

    # Check the condition: {s_1} = C(sigma_1)
    C_result_case2 = C_op(sigma_1_case2)
    
    # Print the "equation" parts
    print(f"The state s_1 is: {s_1_case2}")
    print(f"The set of values sigma_1 is: {sigma_1_case2}")
    print(f"The composed set C(sigma_1) is: {C_result_case2}")
    print(f"The condition {{s_1}} == C(sigma_1) holds: {C_result_case2 == {s_1_case2}}")
    print("\nAs shown, when s_1 != s_0, C(sigma_1) produces multiple states, so the condition fails.")

demonstrate_claim_C()