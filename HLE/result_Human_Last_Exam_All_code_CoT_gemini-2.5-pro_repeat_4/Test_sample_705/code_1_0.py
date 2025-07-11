import itertools
from functools import reduce
from operator import mul

def D_op(states):
    """
    Implements the D operator.
    Decomposes a set of states into a set of their component values.
    """
    decomposed_values = set()
    for state in states:
        decomposed_values.update(state)
    return decomposed_values

def C_op(values, V_sets):
    """
    Implements the C operator.
    Re-composes a set of values into a set of states.
    """
    # For each component k, find the intersection of values with V_k.
    # Per rule 1, if the intersection is empty, all values from V_k are considered possible.
    possible_values_per_component = []
    component_sizes = []
    for V_k in V_sets:
        intersection = values.intersection(V_k)
        if not intersection:
            # Rule 1: D intersects Vk is empty
            possible_values_per_component.append(V_k)
            component_sizes.append(len(V_k))
        else:
            # Rule 2 & 3
            possible_values_per_component.append(intersection)
            component_sizes.append(len(intersection))
    
    # The set of states is the Cartesian product of the possible values for each component.
    recomposed_states = set(itertools.product(*possible_values_per_component))
    
    # Calculate the size of the resulting state set
    size = reduce(mul, component_sizes, 1)

    print(f"Size of C(D) is the product of component possibilities: ", end="")
    print(" * ".join(map(str, component_sizes)), f"= {size}")
    
    return recomposed_states, size

def demonstrate_exponential_blowup():
    """
    Demonstrates that relaxed simulation can require exponentially larger memory.
    """
    # 1. Setup the system
    n = 10  # Number of components
    print(f"Setting up a system with n = {n} components.")
    
    # Each V_k has two values, e.g., V_1 = {'v_1_0', 'v_1_1'}
    V_sets = [set({f'v_{k}_0', f'v_{k}_1'}) for k in range(n)]
    
    # The state space S is the Cartesian product of V_k's. Size is 2^n.
    
    # Define a simulator function f that flips each component value
    # e.g., f( ('v_1_0', 'v_2_0', ...) ) -> ('v_1_1', 'v_2_1', ...)
    def f(state):
        next_state = []
        for k, value in enumerate(state):
            if value == f'v_{k}_0':
                next_state.append(f'v_{k}_1')
            else:
                next_state.append(f'v_{k}_0')
        return tuple(next_state)

    # Choose an initial state s_0
    s_0 = tuple(f'v_{k}_0' for k in range(n))
    print(f"\nOrdinary simulation would track one state at a time, of size {n}.")
    print("Let's trace the relaxed simulation.")

    # 2. Relaxed Simulation: Step 0
    print("\n--- Step 0 ---")
    sigma_0 = D_op({s_0})
    print(f"Initial value set sigma_0 has size: {len(sigma_0)}")
    
    # Compute C(sigma_0) to see the size of the state set we'd simulate from
    print("Calculating the intermediate state set C(sigma_0):")
    C_sigma_0, size_C_sigma_0 = C_op(sigma_0, V_sets)
    assert size_C_sigma_0 == 1
    assert C_sigma_0 == {s_0}

    # 3. Relaxed Simulation: Step 1
    print("\n--- Step 1 ---")
    # First, get the next state from the ordinary simulation to find the new values
    s_1 = f(s_0)
    
    # Update sigma_0 to sigma_1
    new_values = D_op({s_1})
    sigma_1 = sigma_0.union(new_values)
    print(f"After one step, value set sigma_1 has size: {len(sigma_1)}")

    # Now, compute C(sigma_1). This is where the blowup happens.
    # sigma_1 now contains {'v_k_0', 'v_k_1'} for all k.
    print("Calculating the intermediate state set C(sigma_1) for the next step:")
    C_sigma_1, size_C_sigma_1 = C_op(sigma_1, V_sets)
    
    print(f"\nResult: The memory required for the intermediate state set C(sigma_1) is for {size_C_sigma_1} states.")
    print(f"This is exponential in n ({2}^{n}), demonstrating Claim A.")
    print(f"An ordinary simulation step would only handle a single state.")

if __name__ == '__main__':
    demonstrate_exponential_blowup()