import numpy as np
from scipy.stats import linregress
import math

def is_in_cube(pos, n):
    """Checks if a position vector is inside the discrete cube [0, 2n]^3."""
    return np.all((pos >= 0) & (pos <= 2 * n))

def simulate_one_trial(n, outer_box_factor=10):
    """
    Simulates one trial to see if a walk escapes C_n forever.
    Returns True for escape, False for return, None for timeout.
    """
    # Part 1: Walk from X_0 = (n,0,0) until first exit from C_n
    pos = np.array([n, 0, 0])
    
    # max_steps is a safeguard. The walk is transient, so it will exit C_n.
    # The expected time to exit from a point at distance k from the boundary is ~k^2.
    # Here, we start on the boundary, so exit is fast.
    for _ in range(100000): # A generous number of steps
        axis = np.random.randint(3)
        direction = np.random.choice([-1, 1])
        pos[axis] += direction
        
        if not is_in_cube(pos, n):
            # Exited the cube C_n. `pos` is now the first point outside.
            break
    else:
        # Should not be reached in a correct simulation
        return None

    # Part 2: From the exit point, does it return to C_n or escape to "infinity"?
    # "Infinity" is modeled as the boundary of a much larger concentric cube.
    center_of_Cn = np.array([n, n, n])
    infinity_boundary_dist = n * outer_box_factor
    
    # Expected time to travel distance R is ~R^2. This is a safeguard.
    for _ in range(10 * (n * outer_box_factor)**2):
        axis = np.random.randint(3)
        direction = np.random.choice([-1, 1])
        pos[axis] += direction

        if is_in_cube(pos, n):
            return False  # The walk returned to C_n.
        
        # Check if the walk reached our "infinity" boundary
        if np.max(np.abs(pos - center_of_Cn)) > infinity_boundary_dist:
            return True # The walk escaped to infinity.
    
    return None # Timed out, trial is inconclusive.

def main():
    """
    Runs Monte Carlo simulations to estimate p_n and computes the limit
    by fitting a line to the log-log data.
    """
    # We choose n values that are not too small, but also not too large
    # to keep simulation times manageable.
    n_values = [5, 8, 12, 18]
    # Number of trials for each n to get a statistically significant p_n
    num_trials = 30000

    log_n_values = []
    log_inv_p_values = []

    print("Running Monte Carlo simulations to estimate p_n for various n...")
    print("n\t p_n (estimated)\t ln(n)\t\t ln(1/p_n)")
    print("-" * 60)

    for n in n_values:
        successes = 0
        valid_trials = 0
        for _ in range(num_trials):
            result = simulate_one_trial(n)
            if result is not None:
                valid_trials += 1
                if result:
                    successes += 1
        
        if valid_trials > 0:
            p_n = successes / valid_trials
            if p_n > 0:
                log_n = math.log(n)
                log_inv_p = math.log(1 / p_n)
                log_n_values.append(log_n)
                log_inv_p_values.append(log_inv_p)
                print(f"{n}\t {p_n:.6f}\t\t {log_n:.4f}\t\t {log_inv_p:.4f}")
            else:
                print(f"{n}\t {p_n:.6f}\t\t {math.log(n):.4f}\t\t inf (p_n is zero)")
        else:
             print(f"n={n}: No valid trials completed. Consider increasing max_steps.")

    if len(log_n_values) < 2:
        print("\nNot enough data points to perform linear regression.")
        return

    # Using scipy.stats.linregress to fit a line to the log-log data.
    # We are fitting the model: log(1/p_n) = slope * log(n) + intercept
    slope, intercept, r_value, _, _ = linregress(log_n_values, log_inv_p_values)

    print("\n" + "-" * 60)
    print("To find the limit, we fit a line to the points (ln(n), ln(1/p_n)).")
    print("The theoretical relationship is ln(1/p_n) ≈ α * ln(n) + C, where α is the desired limit.")
    print("The slope of the fitted line gives our estimate for α.")
    print("\n--- Linear Regression Result ---")
    print(f"The final fitted equation is: ln(1/p_n) = {slope} * ln(n) + {intercept}")
    print(f"R-squared: {r_value**2}")
    print("--------------------------------\n")
    print(f"The estimated limit is the slope of the line: {slope:.4f}")

if __name__ == '__main__':
    main()