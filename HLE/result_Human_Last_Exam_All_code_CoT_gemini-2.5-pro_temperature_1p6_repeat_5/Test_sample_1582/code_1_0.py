import random

def demonstrate_non_positive_recurrence():
    """
    Simulates a biased random walk to demonstrate non-positive recurrence
    and prints the results.
    """
    # 1. Define the Markov Chain parameters
    p = 0.6  # Probability of moving right (i -> i+1)
    q = 1 - p  # Probability of moving left (i -> i-1)
    num_steps = 1_000_000
    initial_state = 5
    
    # Check the condition from the problem
    # f(x) = x, A = {0}
    # For x > 0, the drift is E[f(X_1)|X_0=x] - f(x) = 2p - 1
    drift = 2 * p - 1

    print("--- Demonstrating Non-Positive Recurrence ---")
    print(f"Using a biased random walk on {{0, 1, 2, ...}} with p={p}.\n")
    
    print("The Lyapunov function is f(x)=x. The finite set is A={0}.")
    print("Checking the condition: E[f(X_1)|X_0=x] - f(x) >= 0 for x > 0.")
    print("The final equation for the drift is:")
    print(f"({p} * (x+1) + {q} * (x-1)) - x = {drift:.2f}")
    if drift >= 0:
        print(f"Since {drift:.2f} >= 0, the condition is satisfied.\n")
    else:
        print(f"Since {drift:.2f} < 0, the condition is NOT satisfied.\n")

    # 2. Simulate the chain
    current_state = initial_state
    path = [current_state]
    state_counts = {}

    for _ in range(num_steps):
        # Update counts
        state_counts[current_state] = state_counts.get(current_state, 0) + 1
        
        # Transition
        if current_state == 0:
            current_state = 1
        else:
            if random.random() < p:
                current_state += 1
            else:
                current_state -= 1
        path.append(current_state)

    # 3. Analyze and print results
    print(f"Simulation run for {num_steps} steps.")
    print(f"Initial state: {initial_state}, Final state: {path[-1]}")
    
    # A positive recurrent chain would have a stationary distribution, meaning it
    # would spend a predictable, non-zero fraction of time in each state.
    # We expect this chain to drift to infinity and spend very little time
    # in the initial states.
    print("\nEmpirical fraction of time spent in states 0-10:")
    for i in range(11):
        count = state_counts.get(i, 0)
        fraction = count / num_steps
        print(f"State {i:2d}: {fraction:.6f} ({count} visits)")

    total_low_state_visits = sum(state_counts.get(i, 0) for i in range(11))
    fraction_low = total_low_state_visits / num_steps
    print(f"\nFraction of time spent in states 0-10: {fraction_low:.4f}")
    print("This low fraction suggests the chain does not repeatedly return to the origin,")
    print("which is inconsistent with positive recurrence.")

demonstrate_non_positive_recurrence()