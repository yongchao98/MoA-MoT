import random

def simulate_random_walk(p, n_steps, start_pos=0):
    """Simulates a random walk on Z_+ with reflection at 0."""
    q = 1.0 - p
    position = start_pos
    path = [position]
    visits_to_A = 1 if start_pos == 0 else 0
    
    for _ in range(n_steps):
        if position == 0:
            position = 1
        else:
            if random.random() < p:
                position += 1
            else:
                position -= 1
        path.append(position)
        if position == 0:
            visits_to_A += 1
            
    return path, visits_to_A

def analyze_case(p, n_steps, n_runs):
    """Runs multiple simulations for a given p and prints analysis."""
    q = 1.0 - p
    
    print(f"\n--- Analysis for p = {p} ---")
    print("The finite set is A = {0} and the function is f(x) = x.")
    
    # We check the condition sum(p(x,y)f(y)) - f(x) >= 0 for x not in A
    # which simplifies to the equation p - q >= 0.
    drift_val = p - q
    print(f"Checking the condition for x > 0: p - q = {p:.1f} - {q:.1f} = {drift_val:.2f}")

    condition_met = drift_val >= 0
    print(f"Is the condition fulfilled? {'Yes' if condition_met else 'No'}")
    
    if condition_met:
        print("Theoretical Conclusion: The chain is NOT positive recurrent.")
    else:
        print("Theoretical Conclusion: The chain could be positive recurrent (and it is).")

    # Run simulations to illustrate the behavior
    total_visits = 0
    final_positions = []
    for i in range(n_runs):
        path, visits_to_A = simulate_random_walk(p, n_steps, start_pos=10) # Start away from A
        final_positions.append(path[-1])
        total_visits += visits_to_A

    avg_final_pos = sum(final_positions) / n_runs
    avg_visits_to_A = total_visits / n_runs

    print("\nSimulation Results:")
    print(f"Average final position after {n_steps} steps: {avg_final_pos:.2f}")
    print(f"Average number of visits to A={{0}} during these steps: {avg_visits_to_A:.2f}")
    
    if condition_met:
        print("Observation: The walker tends to drift away (high final position) and visits A infrequently.")
    else:
        print("Observation: The walker tends to return to the vicinity of A and visits it more frequently.")

if __name__ == '__main__':
    N_STEPS = 5000
    N_RUNS = 100

    # Case 1: Positive Recurrent (Condition NOT met)
    analyze_case(p=0.4, n_steps=N_STEPS, n_runs=N_RUNS)

    # Case 2: Null Recurrent (Condition met)
    analyze_case(p=0.5, n_steps=N_STEPS, n_runs=N_RUNS)

    # Case 3: Transient (Condition met)
    analyze_case(p=0.6, n_steps=N_STEPS, n_runs=N_RUNS)
