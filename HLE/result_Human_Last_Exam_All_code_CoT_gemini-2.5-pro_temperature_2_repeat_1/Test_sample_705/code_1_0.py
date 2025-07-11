import itertools
import sys

def demonstrate_claim_A(n, m):
    """
    Demonstrates the exponential memory growth for relaxed simulation.

    Args:
        n (int): The number of components in the state vector (n in the problem).
        m (int): The number of values to assume are in sigma_i for each component.
    """
    print(f"--- Demonstrating Claim A for a system with n={n} components ---")
    
    # In an ordinary simulation, the memory for computation at each step is
    # primarily for storing one state, which is a tuple of n elements.
    # We can represent its memory size as being proportional to n.
    ordinary_sim_memory_units = n
    print(f"Ordinary Simulation: Stores one state of {n} components.")
    print(f"Computational Memory Requirement is proportional to n = {ordinary_sim_memory_units}.")
    print("-" * 20)

    # In the relaxed simulation, we analyze the computation of sigma_{i+1}.
    # This involves building the set C(sigma_i).
    # Let's assume a scenario where sigma_i has accumulated 'm' distinct values
    # for each of the 'n' components.
    print(f"Relaxed Simulation: Assume sigma_i contains m={m} values for each component.")
    
    # We define the component value sets V_k and the set of known values sigma_i
    # For demonstration, we use integers for values.
    # Vs is the list [V_1, V_2, ..., V_n]
    Vs = [set(range(m * k, m * (k + 1))) for k in range(n)]
    
    # sigma_i is the union of all these component values in this scenario.
    sigma_i = set()
    for v_set in Vs:
        sigma_i.update(v_set)

    # Now, we compute C(sigma_i).
    # According to the rules, if sigma_i has values for a component V_k,
    # those are the options for that component in the Cartesian product.
    # The size of C(sigma_i) is the product of the sizes of these option sets.
    
    # Let's build the list of options for the Cartesian product
    list_of_options = []
    for k in range(n):
        options_k = sigma_i.intersection(Vs[k])
        # This will be Vs[k] itself in our setup, which has size m.
        list_of_options.append(options_k)
        
    size_of_C_set = 1
    for options in list_of_options:
        size_of_C_set *= len(options)
    
    # Each state in C(sigma_i) is a tuple of n elements.
    # So total memory to store the set C(sigma_i) is n * |C(sigma_i)|
    relaxed_sim_memory_units = n * size_of_C_set
    
    print(f"The size of the intermediate set C(sigma_i) to be computed is m^n.")
    equation = f"{m}^{n} = {size_of_C_set}"
    print(f"Equation: {equation}")
    print(f"Memory to store this set is proportional to n * m^n = {n} * {size_of_C_set} = {relaxed_sim_memory_units}.")

    if ordinary_sim_memory_units > 0:
        ratio = relaxed_sim_memory_units / ordinary_sim_memory_units
        print(f"\nThe memory requirement for the relaxed simulation computation is {ratio:.0f} times larger")
        print("than for the ordinary simulation in this case.")
    
    print("\nThis demonstrates that memory requirements for the relaxed simulation can be")
    print("exponentially larger than for the ordinary simulation, supporting Claim A.")


# Run the demonstration for a small system where the difference is clear.
demonstrate_claim_A(n=10, m=2)
