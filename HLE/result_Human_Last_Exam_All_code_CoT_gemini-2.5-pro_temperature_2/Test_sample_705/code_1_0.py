import itertools

def demonstrate_exponential_growth(n):
    """
    This function demonstrates the scenario described in claim A.

    We set up a simulation where:
    - There are 'n' components in a state vector.
    - Each component k can take one of two values from a set Vk.
    - The sets Vk are disjoint, e.g., V_k = {(k, 0), (k, 1)}.
    - The simulator function f flips the value within each component.
    - The ordinary simulation from a starting state s0 cycles between just two states.
    - The relaxed simulation, after one step, requires evaluating f on a set of states
      that is exponentially large in 'n'.
    """

    print(f"--- Demonstration for n = {n} ---")

    # 1. Define the sets Vk
    # Vk = {(k, 0), (k, 1)}
    V = {k: set([(k, 0), (k, 1)]) for k in range(1, n + 1)}

    # Define the full value domain D
    D_union = set()
    for k in V:
        D_union.update(V[k])

    # 2. Define the initial state for the ordinary simulation
    s0 = tuple(sorted(list(D_union)))[::2] # -> ((1,0), (2,0), ...)

    # 3. Define the simulator function f
    # f flips the second element of the pair, e.g. (k, 0) -> (k, 1)
    def f(s):
        return tuple([(k, 1 - val) for k, val in s])

    # 4. Perform ordinary simulation
    s1 = f(s0)
    s2 = f(s1)
    ordinary_trajectory = {s0, s1}
    print("Ordinary Simulation:")
    print(f"s0 = {s0}")
    print(f"s1 = f(s0) = {s1}")
    print(f"s2 = f(s1) = {s2}")
    print(f"Number of unique states in ordinary trajectory: {len(ordinary_trajectory)}")

    # 5. Perform relaxed simulation
    # Step 0:
    sigma_0 = {v for v in s0}

    # C(sigma_0) would just be {s0} because |sigma_0 intersect Vk| = 1 for all k
    states_to_sim_0 = {s0}
    print("\nRelaxed Simulation (Step 0):")
    print(f"sigma_0 = D(s0)")
    print(f"Number of states in C(sigma_0): {len(states_to_sim_0)}")

    # Step 1:
    f_outputs_decomposed = set()
    for s in states_to_sim_0:
        f_outputs_decomposed.update(f(s))
    sigma_1 = sigma_0.union(f_outputs_decomposed)

    print("\nRelaxed Simulation (Step 1):")
    print(f"sigma_1 = sigma_0 U D(f(C(sigma_0))) = {sigma_1}")
    print(f"Is sigma_1 equal to the entire domain D? {sigma_1 == D_union}")
    
    # Now, to compute sigma_2, we need to compute C(sigma_1).
    # Since sigma_1 is D_union, sigma_1 intersect Vk = Vk for all k.
    # The size of C(sigma_1) is the product of |sigma_1 intersect Vk| for all k.
    
    num_states_in_C_sigma_1 = 1
    for k in V:
        # For this setup, intersection_size will be 2
        intersection_size = len(sigma_1.intersection(V[k]))
        num_states_in_C_sigma_1 *= intersection_size
        
    print("\nRelaxed Simulation (Next Step Computation):")
    print(f"To compute sigma_2, we need C(sigma_1).")
    print(f"The number of states in C(sigma_1) is {num_states_in_C_sigma_1}, which is 2^{n}.")

    print("\n--- Conclusion ---")
    print("Memory for ordinary simulation is proportional to a small number of states (2).")
    print(f"Memory for relaxed simulation computation is proportional to the number of states in C(sigma_1), which is 2^{n}.")
    print("This is an exponential difference.")
    print("Therefore, Claim A is correct.")


if __name__ == '__main__':
    # We can demonstrate for a small n like 4. The logic holds for any n.
    demonstrate_exponential_growth(4)
    print("\n")
    demonstrate_exponential_growth(10)
