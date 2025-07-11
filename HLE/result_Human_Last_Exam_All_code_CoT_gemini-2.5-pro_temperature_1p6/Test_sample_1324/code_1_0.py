import math

def chocolate_pass_step(chocolates):
    """
    Performs one step of the chocolate passing game for a list of chocolates.
    Returns the new list of chocolates and the equations for their calculation.
    """
    n = len(chocolates)
    next_chocolates = [0] * n
    equations = []
    
    for k in range(n):
        # Person p_k (at index k) receives from person p_{k-1} (at index k-1).
        # We use modulo for the circular arrangement.
        person_k_val = chocolates[k]
        person_k_minus_1_val = chocolates[(k - 1 + n) % n]
        
        # In the problem, p_b receives from p_a, so b's new chocolates
        # are based on c_b and c_a. Here, p_k is p_b and p_{k-1} is p_a.
        total = person_k_val + person_k_minus_1_val
        half = total / 2
        
        # Note: since all c_k^i are even, their sum is even, and 'half' is always an integer.
        
        equation = f"c_{k+1}^{{i+1}} = (c_{k+1}^i + c_{k}^i) / 2 = ({person_k_val} + {person_k_minus_1_val}) / 2 = {half}"

        if half % 2 != 0:
            # If the result is odd, they take an extra chocolate.
            next_chocolates[k] = int(half) + 1
            equation += f" -> {next_chocolates[k]} (value is odd, add 1)"
        else:
            # If the result is even, it remains as is.
            next_chocolates[k] = int(half)
            equation += f" (value is even)"
            
        equations.append(equation)

    return next_chocolates, equations

def analyze_simulation(initial_chocolates, max_steps):
    """
    Runs the simulation for a given initial state and prints the analysis.
    """
    print(f"Let's analyze the system with n={len(initial_chocolates)} people.")
    chocolates = list(initial_chocolates)
    
    for i in range(max_steps):
        h = max(chocolates)
        l = min(chocolates)
        d = h - l
        
        print(f"\n--- Minute i = {i} ---")
        print(f"State: {chocolates}")
        print(f"h^{i} = {h}, l^{i} = {l}, d^{i} = {d}")
        
        if d == 0:
            print("The system has reached a stable state where all values are equal.")
            break
        
        print(f"\nCalculating the state for minute i+1 = {i+1}:")
        new_chocolates, equations = chocolate_pass_step(chocolates)
        for eq in equations:
            print(eq)
        chocolates = new_chocolates
        
    # Print the final state reached
    if d > 0:
      print(f"\n--- After {max_steps} minutes ---")
      h = max(chocolates)
      l = min(chocolates)
      d = h - l
      print(f"State: {chocolates}")
      print(f"h = {h}, l = {l}, d = {d}")

# --- Main Execution ---

# Case 1: A dynamic example that evolves over time.
initial_state_1 = [10, 4, 2, 2]
analyze_simulation(initial_state_1, 8)

# Case 2: An example that starts in a stable (converged) state.
print("\n" + "="*50 + "\n")
initial_state_2 = [8, 8, 8, 8]
analyze_simulation(initial_state_2, 3)
