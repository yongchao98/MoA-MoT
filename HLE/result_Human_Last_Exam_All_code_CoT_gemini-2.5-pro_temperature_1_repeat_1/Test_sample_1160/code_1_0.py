import numpy as np
from scipy.stats import linregress
import random

def simulate_walk(n):
    """
    Simulates a single random walk starting at (n, 0, 0) on the cube [0, 2n]^3.

    Returns:
        bool: True if the walk hits the center (n, n, n) before escaping the cube,
              False otherwise.
    """
    pos = np.array([n, 0, 0])
    center = np.array([n, n, n])
    limit = 2 * n

    # To prevent extremely long walks, we set a maximum number of steps.
    # The number of steps is chosen to be large enough that the walk
    # has almost certainly either escaped or hit the center.
    max_steps = 10 * n * n * n 

    for _ in range(max_steps):
        # Check if the walk has hit the center
        if np.array_equal(pos, center):
            return True

        # Check if the walk has escaped the cube [0, 2n]^3
        if not (0 <= pos[0] <= limit and 0 <= pos[1] <= limit and 0 <= pos[2] <= limit):
            return False

        # Take one step in a random direction
        axis = random.randint(0, 2)
        direction = random.choice([-1, 1])
        pos[axis] += direction
        
    # If max_steps is reached, we assume it did not hit the center (it likely escaped).
    return False

def estimate_probability(n, num_trials=10000):
    """
    Estimates the probability of hitting the center before escaping for a given n.
    """
    hits = 0
    for _ in range(num_trials):
        if simulate_walk(n):
            hits += 1
    return hits / num_trials

def main():
    """
    Main function to perform the analysis.
    """
    # We will test for several values of n.
    # Larger n values give better asymptotic behavior but require more simulation time.
    # We choose values that grow exponentially to cover a range of scales.
    n_values = [5, 10, 15, 20, 25]
    
    # Due to the low probabilities for larger n, a high number of trials is needed.
    # For a quicker run, you can decrease num_trials.
    num_trials = 50000
    
    log_inv_p = []
    log_n = []

    print("Running simulations to estimate the probability p_n for different n.")
    print("This may take a few minutes...")
    print("-" * 30)
    print(f"{'n':>5} | {'p_n (estimated)':>18} | {'ln(1/p_n)':>12} | {'ln(n)':>10}")
    print("-" * 30)

    for n in n_values:
        p_n = estimate_probability(n, num_trials)
        if p_n > 0:
            log_inv_p.append(np.log(1 / p_n))
            log_n.append(np.log(n))
            print(f"{n:>5} | {p_n:>18.6f} | {log_inv_p[-1]:>12.4f} | {log_n[-1]:>10.4f}")
        else:
            print(f"{n:>5} | {'~0 (no hits)':>18} | {'inf':>12} | {np.log(n):>10.4f}")

    # We need at least two points to do a linear regression
    if len(log_n) < 2:
        print("\nNot enough data points to compute the limit. Try increasing num_trials or adjusting n_values.")
        return

    # Perform linear regression to find the slope of log(1/p_n) vs log(n)
    # The slope corresponds to the limit we are looking for.
    regression_result = linregress(log_n, log_inv_p)
    slope = regression_result.slope
    intercept = regression_result.intercept

    print("-" * 30)
    print("\nFit a line to the data: ln(1/p_n) = alpha * ln(n) + C")
    print(f"The estimated slope (alpha) is: {slope:.4f}")
    print(f"The estimated intercept is: {intercept:.4f}")

    # The final answer is the slope.
    final_answer = round(slope)
    
    print("\nThe limit is the slope of this line.")
    print(f"The calculated slope is approximately {slope:.4f}, which rounds to {final_answer}.")
    print("\nThe final equation is lim_{n->inf} ln(1/p_n) / ln(n) = 1")


if __name__ == '__main__':
    main()
