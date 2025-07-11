import random

def demonstrate_non_positive_recurrence():
    """
    Simulates a Markov chain that satisfies the given conditions and
    demonstrates that it is not positive recurrent.
    """
    # Parameters for a birth-death chain on N_0 = {0, 1, ...}
    # Let p_i = p(i, i+1) and q_i = p(i, i-1).
    # We choose p > q to create a drift to infinity.
    p = 0.6  # Probability of i -> i+1 for i > 0
    q = 0.4  # Probability of i -> i-1 for i > 0

    # Simulation parameters
    n_steps = 10000
    start_state = 5
    
    # Define the problem setup for this example
    # A is a finite set. Let A = {0}.
    # f(x) is a function. Let f(x) = x. f(x) -> infinity as x -> infinity.
    print("--- Problem Setup ---")
    print(f"State space: Non-negative integers {{0, 1, 2, ...}}")
    print(f"Finite set A = {{0}}")
    print(f"Function f(x) = x")
    print(f"Transition probabilities for x > 0: p(x, x+1) = {p}, p(x, x-1) = {q}")
    print(f"At state 0, let's assume it always moves to 1 (reflecting barrier).")
    
    # --- Verify the Lyapunov-like condition ---
    # For all x not in A (i.e., x > 0), we need E[f(X_1) | X_0=x] - f(x) >= 0.
    # The equation is: p*f(x+1) + q*f(x-1) - f(x) >= 0
    # Substituting f(x)=x: p*(x+1) + q*(x-1) - x >= 0
    # -> px + p + qx - q - x >= 0
    # -> (p+q)x + p-q - x >= 0
    # -> x + p-q - x >= 0
    # -> p - q >= 0
    print("\n--- Verifying the Condition from the Problem ---")
    print("For x not in A, the condition is: E[f(X_1)|X_0=x] - f(x) >= 0")
    print("With f(x)=x and our probabilities, this simplifies to the equation: p - q >= 0")
    print("Let's plug in the numbers:")
    p_minus_q = p - q
    print(f"{p} - {q} = {p_minus_q:.2f}")
    is_satisfied = p_minus_q >= 0
    print(f"Is the condition satisfied? {is_satisfied}")
    
    if not is_satisfied:
        print("\nWarning: The chosen p and q do not satisfy the condition.")
        return

    # --- Simulation to show behavior ---
    print(f"\n--- Simulating the chain for {n_steps} steps from x_0 = {start_state} ---")
    current_state = start_state
    path = [current_state]
    visits_to_A = 0
    
    for _ in range(n_steps):
        if current_state == 0:
            # We are in A. Let's move out.
            current_state = 1
        else:
            # We are outside A. Standard birth-death step.
            if random.random() < p:
                current_state += 1
            else:
                current_state -= 1
        path.append(current_state)
        if current_state == 0:
            visits_to_A += 1
            
    print("\n--- Simulation Results ---")
    print(f"Final state after {n_steps} steps: {path[-1]}")
    print(f"Number of visits to A={{0}} during the simulation: {visits_to_A}")

    # A positive recurrent chain would return to A frequently.
    # The low number of visits suggests the chain is not positive recurrent.
    if visits_to_A > 0:
        visit_frequency = (n_steps + 1) / visits_to_A
        print(f"Approximate mean time between visits to A: {visit_frequency:.2f} steps.")
    else:
        print("Set A was never revisited after leaving the starting region.")
        
    print("\nThe drift to infinity means the state tends to increase over time, as seen below.")
    print("Path (first 20 and last 20 states):")
    if n_steps <= 40:
        print(path)
    else:
        print(f"Start: {path[:20]}...")
        print(f"End:   ...{path[-20:]}")
        

demonstrate_non_positive_recurrence()
