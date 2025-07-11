import numpy as np
from scipy.stats import linregress

def estimate_fluctuation_scaling():
    """
    Numerically simulates the ODE to find the scaling of fluctuation magnitude R with epsilon.
    """
    # Epsilon values to test. Using a logarithmic scale is best for finding scaling laws.
    epsilons = np.logspace(-2.2, -1.2, 6)
    r_values = []

    print("Running simulations for different epsilon values...")

    for epsilon in epsilons:
        L = 1.0 / epsilon
        N = int(L) - 1
        if N <= 0:
            continue

        num_trials = 200  # Number of random simulations for each epsilon
        num_x_points = 500 # Resolution for the x-domain

        x_grid = np.linspace(0, L, num_x_points, endpoint=False)[1:] # Grid for x, avoiding boundaries

        # Store the result of y2s for each trial
        y2s_trials = np.zeros((num_trials, len(x_grid)))

        for i in range(num_trials):
            # Generate N ordered random points z_i in [0, L]
            z = np.sort(np.random.uniform(0, L, N))

            # Calculate y2s(x) = sum_i G(x, z_i) using vectorization
            # G(x, z) = x(z-L)/L for x < z and z(x-L)/L for x > z
            x_col = x_grid[:, np.newaxis]
            z_row = z[np.newaxis, :]
            
            G = np.where(x_col < z_row, x_col * (z_row - L) / L, z_row * (x_col - L) / L)
            y2s = np.sum(G, axis=1)
            y2s_trials[i, :] = y2s

        # Calculate variance of y2s at each x point across all trials
        var_y2s = np.var(y2s_trials, axis=0)
        
        # Find the maximum variance
        max_var_y2s = np.max(var_y2s)

        # Calculate R = epsilon^2 * sqrt(max_var(y2s))
        R = epsilon**2 * np.sqrt(max_var_y2s)
        r_values.append(R)
        print(f"  epsilon = {epsilon:.4f}, N = {N}, R = {R:.6f}")

    # Perform a log-log linear regression to find the scaling exponent
    log_epsilons = np.log(epsilons)
    log_r_values = np.log(r_values)

    # slope is the scaling exponent, intercept is log(C)
    slope, intercept, r_value, p_value, std_err = linregress(log_epsilons, log_r_values)
    
    C = np.exp(intercept)

    print("\n--- Results ---")
    print(f"Log-log regression slope (scaling exponent p): {slope:.4f}")
    print(f"Log-log regression coefficient (C): {C:.4f}")
    print("\nThe estimated relationship is of the form R = C * epsilon^p.")
    print("\nFinal Equation:")
    print(f"R = {C:.4f} * epsilon^{slope:.4f}")

estimate_fluctuation_scaling()
<<<0.5>>>