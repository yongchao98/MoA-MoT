import random

def simulate_biased_random_walk(p, start_state, num_steps):
    """
    Simulates a biased random walk on non-negative integers.

    Args:
        p (float): Probability of moving from x to x+1 (for x > 0).
        start_state (int): The starting state of the walker.
        num_steps (int): The total number of steps to simulate.

    Returns:
        None. Prints the state at regular intervals.
    """
    if not 0.5 <= p <= 1.0:
        print("Warning: p should be >= 0.5 to satisfy the problem's conditions.")

    current_state = start_state
    visits_to_A = 0
    print(f"Starting simulation with p = {p}, start_state = {start_state}, for {num_steps} steps.")
    print("-" * 30)

    print(f"Time 0: State = {current_state}")
    
    report_interval = num_steps // 10

    for step in range(1, num_steps + 1):
        if current_state == 0:
            # From state 0 (our set 'A'), always move to 1
            current_state = 1
        else:
            # For x > 0, move to x+1 with probability p, and x-1 with 1-p
            if random.random() < p:
                current_state += 1
            else:
                current_state -= 1
        
        # Count visits to the set A = {0}
        if current_state == 0:
            visits_to_A += 1

        # Report status at intervals
        if step % report_interval == 0:
            print(f"Time {step}: State = {current_state}")

    print("-" * 30)
    print(f"Simulation finished.")
    print(f"Final state: {current_state}")
    print(f"Total number of visits to A={{0}}: {visits_to_A} in {num_steps} steps.")
    print("\nAs the simulation shows, the state tends to increase over time, drifting away from the origin.")
    print("This behavior is characteristic of a chain that is not positive recurrent.")

if __name__ == '__main__':
    # Parameters for the simulation
    PROB_UP = 0.6  # p > 0.5, so the chain is transient (not positive recurrent)
    START_STATE = 5
    NUM_STEPS = 50000

    simulate_biased_random_walk(PROB_UP, START_STATE, NUM_STEPS)
