import random
import numpy as np

def simulate_walk(n, max_steps):
    """
    Simulates a single random walk and determines if it escapes or returns.
    
    Args:
        n: The parameter defining the cube size.
        max_steps: Maximum number of steps to prevent excessively long walks.

    Returns:
        True if the walk escapes the cube before returning to the start.
        False if the walk returns to the start first.
        None if max_steps is reached.
    """
    start_pos = (n, 0, 0)
    pos = list(start_pos)
    
    # The first step is away from the starting position
    
    # Randomly choose a dimension to move in (x, y, or z)
    dim = random.randint(0, 2)
    # Randomly choose a direction (+1 or -1)
    move = random.choice([-1, 1])
    pos[dim] += move

    for _ in range(max_steps):
        # Check for returning to the starting point
        if tuple(pos) == start_pos:
            return False  # Returned
            
        # Check for escaping the cube
        if not (0 <= pos[0] <= 2 * n and 0 <= pos[1] <= 2 * n and 0 <= pos[2] <= 2 * n):
            return True  # Escaped

        # Perform the next step of the walk
        dim = random.randint(0, 2)
        move = random.choice([-1, 1])
        pos[dim] += move
        
    return None # Timed out

def estimate_pn(n, num_trials, max_steps_factor=20):
    """
    Estimates the escape probability p_n for a given n.
    
    Args:
        n: The parameter defining the cube size.
        num_trials: The number of walks to simulate for the estimation.
        max_steps_factor: Factor to determine max_steps based on n^2.

    Returns:
        The estimated escape probability p_n.
    """
    successes = 0
    timeouts = 0
    max_steps = max_steps_factor * (2 * n)**2
    for _ in range(num_trials):
        result = simulate_walk(n, max_steps)
        if result is True:
            successes += 1
        elif result is None:
            timeouts += 1

    valid_trials = num_trials - timeouts
    if valid_trials > 0:
        return successes / valid_trials
    else:
        return 0.0

def main():
    """
    Main function to run the analysis.
    """
    # Values of n to test
    ns = np.array([5, 10, 15, 20, 25])
    num_trials = 200000  # Increase for better accuracy, but longer runtime
    
    pns = []
    print("Running simulations...")
    for n in ns:
        print(f"  Estimating p_n for n = {n}...")
        p = estimate_pn(n, num_trials)
        pns.append(p)
    pns = np.array(pns)
    
    print("\n--- Results ---")
    
    # Filter out cases where probability is zero to avoid log errors
    valid_indices = pns > 0
    if not np.any(valid_indices):
        print("All simulations resulted in p_n=0. Cannot perform analysis.")
        print("This could be due to too few trials or walks that are too long.")
        return
        
    final_ns = ns[valid_indices]
    final_pns = pns[valid_indices]

    # Calculate log values
    log_ns = np.log(final_ns)
    log_1_over_pns = np.log(1.0 / final_pns)
    
    # Perform linear regression
    slope, intercept = np.polyfit(log_ns, log_1_over_pns, 1)

    print("The final equation is of the form: log(1/p_n) = a * log(n) + b")
    print("Values for log(n):")
    for val in log_ns:
        print(f"{val:.4f}")
    
    print("\nValues for log(1/p_n):")
    for val in log_1_over_pns:
        print(f"{val:.4f}")

    print(f"\nEstimated slope (a): {slope:.4f}")
    print(f"Estimated intercept (b): {intercept:.4f}")
    
    print("\nThe theoretical limit is 1. The simulation result is an approximation.")

if __name__ == "__main__":
    main()
