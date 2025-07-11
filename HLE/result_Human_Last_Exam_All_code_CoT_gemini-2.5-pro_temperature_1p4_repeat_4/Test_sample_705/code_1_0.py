import math

def demonstrate_exponential_blowup():
    """
    This function demonstrates that relaxed simulation can require exponentially
    larger computational space than ordinary simulation.
    """
    n = 10  # Number of components in the state vector

    # 1. Define the disjoint sets V_k
    # V_k = {10*k, 10*k + 1}. e.g., V_1={10,11}, V_2={20,21}
    V_sets = [set([10 * k, 10 * k + 1]) for k in range(1, n + 1)]

    # 2. Define an initial state s_0
    # s_0 = (10, 20, ..., 100)
    s0 = tuple(10 * k for k in range(1, n + 1))
    
    # 3. Define a simulator function f(s)
    # This function is designed to introduce new values.
    # It maps s0 to a state with all the other possible values.
    s_other = tuple(10 * k + 1 for k in range(1, n + 1))
    def f(s):
        if s == s0:
            return s_other
        # Behavior for other states doesn't matter for this example
        return s

    print(f"System parameters: n = {n}")
    print(f"Initial state s_0 = {s0}\n")

    # --- Ordinary Simulation ---
    # At each step, we only need to compute f for ONE state.
    print("--- Ordinary Simulation Analysis ---")
    s1 = f(s0)
    print("To compute s_1 = f(s_0), we process 1 state.")
    # To compute s_2 = f(s_1), we would again process 1 state.
    # The computational memory requirement is proportional to a single state.
    
    # --- Relaxed Simulation ---
    print("\n--- Relaxed Simulation Analysis ---")
    # Step 0 -> 1
    # sigma_0 is the set of values in s_0
    sigma_0 = set(s0)
    # C(sigma_0) has only one state, s_0 itself.
    # We compute f(s_0) to get values for the next step.
    s1_relaxed = f(s0)
    
    # sigma_1 = sigma_0 U D(f(s_0))
    sigma_1 = sigma_0.union(set(s1_relaxed))
    
    print(f"sigma_0 = {sigma_0}")
    print(f"s_1 = f(s_0) = {s1_relaxed}")
    print(f"sigma_1 = sigma_0 U D(s_1) = {sigma_1}")
    
    # Step 1 -> 2
    # To compute sigma_2, we must first compute C(sigma_1).
    # The size of C(sigma_1) determines the number of f() evaluations needed.
    # Size of C(sigma_1) = product over k of |sigma_1 intersect V_k|
    
    # Calculate the size of the set of states to process
    component_sizes = []
    for v_k in V_sets:
        intersection_size = len(sigma_1.intersection(v_k))
        component_sizes.append(intersection_size)
    
    num_states_to_process = math.prod(component_sizes)

    print("\nTo compute sigma_2, we need to process all states in C(sigma_1).")
    print("The number of states to process is the size of C(sigma_1).")

    # Build the equation string
    equation_str = " * ".join(map(str, component_sizes))
    
    print(f"Size of C(sigma_1) = {equation_str} = {num_states_to_process}")

    print("\nConclusion:")
    print("Ordinary simulation step requires processing 1 state.")
    print(f"Relaxed simulation step requires processing {num_states_to_process} states.")
    print(f"This demonstrates an exponential ({2}^{n}) growth in computational (and memory) requirements for the relaxed simulation.")

demonstrate_exponential_blowup()