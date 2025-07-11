def simulate_chocolate_game(initial_chocolates, steps):
    """
    Simulates the chocolate passing game for a given number of steps.

    This function implements the rules described and prints the state at each step,
    including the explicit calculation based on Equation (1).

    Args:
        initial_chocolates (list): A list of even integers for the initial state.
        steps (int): The number of minutes (iterations) to simulate.
    """
    chocolates = list(initial_chocolates)
    n = len(chocolates)

    if not all(c % 2 == 0 for c in chocolates):
        print("Error: All initial chocolate counts must be even.")
        return

    print("--- Chocolate Game Simulation ---")
    print(f"Number of people (n): {n}\n")
    
    # Print initial state
    print(f"Initial state (i=0):")
    print(f"  Chocolates c^0: {chocolates}")
    h = max(chocolates)
    l = min(chocolates)
    d = h - l
    print(f"  h^0 = {h}, l^0 = {l}, d^0 = {d}\n")

    for i in range(1, steps + 1):
        prev_chocolates = list(chocolates)
        next_chocolates = [0] * n

        print(f"State after minute {i} (i={i}):")
        print(f"  Calculations for c^{i}:")
        
        # Each person k receives from person (k-1)
        for k in range(n):
            prev_person_idx = (k - 1 + n) % n
            c_a = prev_chocolates[prev_person_idx] # Person giving chocolates
            c_b = prev_chocolates[k]           # Person receiving chocolates
            
            # Temporary value is the average
            temp_val = (c_a + c_b) // 2

            # Apply Equation (1)
            final_val = 0
            is_odd = (temp_val % 2 != 0)
            if is_odd:
                final_val = temp_val + 1
            else:
                final_val = temp_val

            next_chocolates[k] = final_val

            # Print the calculation for this person
            # The person index is k, but for readability p_1, p_2, etc.
            print(f"    For p_{k+1}: c_{k+1}^{i} = ", end="")
            if not is_odd:
                print(f"({c_a} + {c_b})/2 = {temp_val} (even)")
            else:
                print(f"({c_a} + {c_b})/2 + 1 = {temp_val} + 1 = {final_val} (since {temp_val} is odd)")
        
        chocolates = next_chocolates
        
        # Print summary for the new state
        print(f"  Chocolates c^{i}: {chocolates}")
        h = max(chocolates)
        l = min(chocolates)
        d = h - l
        print(f"  h^{i} = {h}, l^{i} = {l}, d^{i} = {d}\n")
        
        if d == 0:
            print("System reached equilibrium. Halting simulation.")
            break

# Let's run a simulation with an example that shows the dynamics
n_people = 4
initial_config = [12, 12, 8, 8] 
simulation_steps = 5

simulate_chocolate_game(initial_config, simulation_steps)
