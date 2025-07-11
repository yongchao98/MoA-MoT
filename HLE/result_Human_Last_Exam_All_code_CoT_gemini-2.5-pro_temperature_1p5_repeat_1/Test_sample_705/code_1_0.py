def explain_and_calculate(n, m):
    """
    Analyzes claim A by constructing a specific scenario and calculating the
    computational load for both simulation types.

    Args:
        n (int): The number of disjoint sets (dimensionality of the state space).
        m (int): The number of elements in each set V_k.
    """

    print("--- Analysis of Claim A ---")
    print("We demonstrate that for a specific setup, the relaxed simulation can require an exponentially larger computational load than the ordinary simulation.")
    print("We interpret 'memory space for computation' or 'computational load' as the number of states to simulate at a given step.")
    
    print(f"\n[Setup]")
    print(f"Let n = {n} (number of components).")
    print(f"Let each V_k be a disjoint set of size |V_k| = {m}.")
    print(f"For this example, let's say V_k = {{v_k0, v_k1, ...}}.")

    print("\n[Ordinary Simulation]")
    print("At each step i, the simulation computes s_{i+1} = f(s_i).")
    print("This process always involves a single state evaluation.")
    print("Computational Load: 1 state")

    print("\n[Relaxed Simulation]")
    print("--- Step 0 ---")
    print("Let the initial state be s_0 = (v_10, v_20, ..., v_n0).")
    print("The initial value set is sigma_0 = D({s_0}) = {v_10, v_20, ..., v_n0}.")
    print("To compute sigma_1, we need to evaluate f(s) for all s in C(sigma_0).")
    print("For each component k, the intersection |sigma_0 intersect V_k| = 1.")
    print(f"Thus, the number of states to simulate is |C(sigma_0)| = {' * '.join(['1']*n)} = 1.")
    print("At this step, the computational load is the same as the ordinary simulation.")

    print("\n--- Step 1 ---")
    print("Let's choose a function f such that the next state s_1 = f(s_0) introduces a new value for each component.")
    print("For example, s_1 = (v_11, v_21, ..., v_n1).")
    print("The next value set is sigma_1 = sigma_0 union D({s_1}).")
    print("Now, for each component k, the intersection |sigma_1 intersect V_k| = |{v_k0, v_k1}| = 2.")
    print("To compute sigma_2, we must evaluate f(s) for all s in C(sigma_1).")
    print("The number of states to simulate is |C(sigma_1)|, calculated as the product of the sizes of these intersections.")
    
    # Calculate the size of C(sigma_1)
    sizes_of_intersections = [m] * n
    total_states_to_simulate = 1
    for size in sizes_of_intersections:
        total_states_to_simulate *= size
        
    equation_parts = [str(s) for s in sizes_of_intersections]
    equation_string = " * ".join(equation_parts)

    print(f"|C(sigma_1)| = {equation_string} = {total_states_to_simulate}")

    print("\n[Conclusion]")
    print(f"At step 1, the relaxed simulation must consider {total_states_to_simulate} states, while the ordinary simulation still only considers 1.")
    print(f"The ratio of computational load is {total_states_to_simulate}, which equals {m}^{n}.")
    print(f"This number grows exponentially with n, confirming Claim A.")


# Run the demonstration for a clear example.
# n=10 and m=2 will show an exponential explosion.
explain_and_calculate(n=10, m=2)