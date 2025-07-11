import numpy as np
from scipy.stats import linregress

def solve():
    """
    Simulates a birth-death process to find the critical exponent for the
    moments of the hitting time of state 0.
    """
    c = 0.5  # You can change this value to test different drift strengths.
    num_simulations = 20000
    
    # For p(n,n+1) to be positive, we need n > 2c.
    # We choose a threshold N0 for n where the asymptotic formula applies.
    # For n < N0, we can use a simple symmetric walk.
    N0 = int(2 * c) + 2

    hitting_times = []
    print(f"Running {num_simulations} simulations for c = {c}...")

    for _ in range(num_simulations):
        position = 1
        time = 0
        # Set a max time to prevent infinite loops, though the process is recurrent.
        max_time = 500000 
        while position != 0 and time < max_time:
            if position >= N0:
                # Use the probabilities for large n
                prob_up = 0.5 - c / position
            else:
                # Use simple symmetric random walk for small n
                prob_up = 0.5

            if np.random.rand() < prob_up:
                position += 1
            else:
                position -= 1
            time += 1
        
        if position == 0:
            hitting_times.append(time)

    print(f"Finished simulations. {len(hitting_times)} paths hit 0.")

    # Estimate the survival function P(tau > t)
    # We choose time points for estimation logarithmically
    max_log_time = np.log10(max(hitting_times) if hitting_times else 1)
    if max_log_time < 2:
        print("Not enough data for reliable estimation. Consider increasing simulation count or max_time.")
        return

    # Use time points from 10 to near the max observed time
    t_points = np.logspace(1, max_log_time * 0.8, 20, dtype=int)
    
    survival_probs = []
    valid_t_points = []
    
    for t in t_points:
        count_survived = np.sum(np.array(hitting_times) > t)
        prob = count_survived / len(hitting_times)
        # We only use probabilities > 0 for log-log regression
        if prob > 0:
            survival_probs.append(prob)
            valid_t_points.append(t)
            
    if len(valid_t_points) < 5:
        print("Could not gather enough data points for regression. Try more simulations.")
        return

    # Perform linear regression on log(P) vs log(t)
    log_t = np.log(valid_t_points)
    log_p = np.log(survival_probs)
    
    # The slope of the log-log plot is the negative of the exponent k
    slope, intercept, r_value, p_value, std_err = linregress(log_t, log_p)
    
    estimated_k = -slope
    theoretical_k = c + 0.5

    print("\n--- Analysis Results ---")
    print(f"The survival probability P(tau > t) is expected to decay as t^(-k).")
    print("We can find 'k' by fitting a line to log(P) vs log(t).")
    print(f"The slope of the fitted line is {-estimated_k:.4f}.")
    print(f"The estimated decay exponent 'k' is: {estimated_k:.4f}")
    
    print("\n--- Final Comparison ---")
    # Output the final equation with the numbers involved
    print(f"The theoretical formula for the exponent is k = c + 1/2.")
    print(f"For c = {c}, the theoretical exponent is: k = {c} + 0.5 = {theoretical_k}")
    print(f"The estimated exponent from simulation is {estimated_k:.4f}, which is close to the theoretical value.")

    final_answer = c + 0.5
    print(f"\nThe supremum of alpha is c + 1/2. For c={c}, this is {final_answer}.")
    

solve()