import random

def simulate_chain(steps, p_right):
    """
    Simulates a biased random walk on non-negative integers.

    Args:
        steps (int): The number of steps to simulate.
        p_right (float): The probability of moving to the right (x -> x+1).

    Returns:
        float: The proportion of time spent at state 0.
    """
    position = 0
    time_at_origin = 0
    p_left = 1 - p_right

    for i in range(steps):
        if position == 0:
            # At origin, must move to 1
            position = 1
            # We count the time spent at the origin before the move
            time_at_origin += 1
        else:
            # Biased random walk for x > 0
            if random.random() < p_right:
                position += 1
            else:
                position -= 1
    
    return time_at_origin / steps

def main():
    """
    Main function to run the simulation and print the conclusion.
    """
    # --- Setup ---
    # Number of simulation steps
    num_steps = 1000000
    # Probability of moving right. p > 0.5 ensures drift away from the origin.
    p = 0.6 

    # --- Verification of conditions ---
    # Markov Chain on Sigma = {0, 1, 2, ...}
    # Transition p(x, x+1) = p, p(x, x-1) = 1-p for x > 0. p(0,1) = 1.
    # Let A = {0}, which is finite.
    # Let f(x) = x. f is non-negative and f(x) -> infinity as x -> infinity.
    # For x not in A (i.e., x >= 1):
    # E[f(X_1) | X_0=x] - f(x) = [p * (x+1) + (1-p) * (x-1)] - x
    #                        = (px + p + x - 1 - px + p) - x
    #                        = (x + 2p - 1) - x
    #                        = 2p - 1
    # Since p = 0.6, 2*p - 1 = 1.2 - 1 = 0.2, which is >= 0.
    # All conditions of the problem are satisfied.

    print(f"Simulating a Markov chain for {num_steps} steps.")
    print(f"The chain has a drift to the right (p={p}).")
    print("This chain satisfies the conditions from the problem statement.")
    print("We will calculate the proportion of time the chain spends at state 0.")
    print("For a positive recurrent chain, this should converge to a positive number.")
    print("For a transient or null recurrent chain, this will converge to 0.")
    
    proportion = simulate_chain(num_steps, p)

    print("\n--- Simulation Result ---")
    print(f"Proportion of time spent at state 0: {proportion:.6f}")

    if proportion < 1e-4:
        print("\nConclusion: The very small proportion suggests the chain does not return to the origin often enough.")
        print("This is consistent with a transient or null recurrent chain, not a positive recurrent one.")
    else:
        print("\nConclusion: The result is inconclusive from this simulation alone, but theory predicts non-positive recurrence.")

main()
