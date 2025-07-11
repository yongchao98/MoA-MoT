import numpy as np
from scipy.stats import linregress
import warnings

def estimate_fluctuation_scaling():
    """
    Main function to estimate the scaling of R(epsilon).
    It solves the ODE for many random realizations of the delta function sources
    and calculates the variance of the solution. By doing this for several
    epsilon values, it can determine the scaling law R(epsilon).
    """
    # Epsilon values to test. They should be small and varied.
    epsilons = np.array([0.1, 0.08, 0.05, 0.04, 0.02, 0.01])
    # Number of Monte Carlo trials for each epsilon to get good statistics
    num_trials = 200
    # Number of grid points for the spatial variable x
    num_x_points = 2000
    
    R_values = []

    print("Running simulations for different epsilon values...")
    # Loop over each chosen epsilon
    for eps in epsilons:
        print(f"  epsilon = {eps:.3f}")
        L = 1.0 / eps
        N = int(L) - 1
        x_grid = np.linspace(0, L, num_x_points)
        
        # Array to store the results of all trials for this epsilon
        y_trials = np.zeros((num_trials, num_x_points))

        # Loop for Monte Carlo trials
        for i in range(num_trials):
            # As per the problem, z_i are ordered values from a uniform distribution.
            # This corresponds to the order statistics of i.i.d. uniform variables.
            z_i_random = np.random.uniform(0, L, size=N)
            z_i_ordered = np.sort(z_i_random)
            
            # Analytically solve the ODE for this specific realization of z_i
            y_trials[i, :] = solve_for_y_realization(x_grid, eps, z_i_ordered)

        # The variance of the fluctuation is Var[y(x) - y0(x)].
        # Since y0(x) is the (deterministic) mean of y(x), this is equal to Var[y(x)].
        var_y = np.var(y_trials, axis=0)
        
        # R^2 is the maximum of this variance over x
        max_var = np.max(var_y)
        
        # R is the square root of the max variance
        R = np.sqrt(max_var)
        R_values.append(R)

    # Fit a power law R = C * epsilon^p by performing linear regression on the logs
    log_eps = np.log(epsilons)
    log_R = np.log(R_values)
    
    # Use scipy's linregress to find the slope (p) and intercept (log C)
    p, log_C, _, _, _ = linregress(log_eps, log_R)
    C = np.exp(log_C)
    
    print("\n--- Results ---")
    print(f"The analysis has been performed over {len(epsilons)} values of epsilon with {num_trials} trials each.")
    print(f"The simulation data is fitted to the power law R(epsilon) = C * epsilon^p.")
    print(f"The obtained scaling exponent is p = {p:.4f}")
    print(f"The obtained prefactor is C = {C:.4f}")
    print("\nThe estimated relationship between the fluctuation magnitude R and epsilon is:")
    # Final formatted output of the equation
    print(f"R(epsilon) = {C:.4f} * epsilon^{p:.4f}")
    return p

def solve_for_y_realization(x_grid, epsilon, z):
    """
    Solves the ODE for a single given realization of source locations z_i.
    The analytical solution is:
    y(x) = e^(ex) * (1 + Integral_0^x [C + e^2 N(t)]e^(-et) dt)
    where N(t) is the number of z_i's less than t, and C is a constant
    determined by the boundary condition y(L)=0.
    """
    L = 1.0 / epsilon
    N = z.size
    
    # The counting function N(t) is a step function. The integral involving it
    # can be computed exactly by summing over the intervals between z_i's.
    z_full = np.concatenate(([0], z, [L]))
    
    # First, calculate Integral_0^L N(t)exp(-et) dt to find C
    integral_N_L = 0.0
    for i in range(N + 1):
        # In the interval [z_i, z_{i+1}), N(t) is constant and equal to i
        # The integral of i*exp(-et) is -(i/e)*exp(-et)
        val1 = np.exp(-epsilon * z_full[i])
        val2 = np.exp(-epsilon * z_full[i+1])
        integral_N_L += (i / epsilon) * (val1 - val2)

    # The constant C is y'(0) - epsilon*y(0). We find it from y(L)=0.
    # The condition gives: 1 + C * integral(e^-et) + e^2 * integral(N(t)e^-et) = 0
    C = -epsilon * (1 + epsilon**2 * integral_N_L) / (1 - np.exp(-1.0))

    # Now, calculate y(x) for each x in the grid.
    # This involves calculating Integral_0^x for each x.
    integral_C_term = C * (1 - np.exp(-epsilon * x_grid)) / epsilon
    
    integral_N_term = np.zeros_like(x_grid, dtype=float)
    # This integral is calculated piecewise.
    for i in range(N + 1):
        count = i
        lower_bound = z_full[i]
        upper_bound = z_full[i+1]
        
        # Add contribution for all x values beyond the current interval
        idx = x_grid > upper_bound
        integral_N_term[idx] += (count / epsilon) * (np.exp(-epsilon * lower_bound) - np.exp(-epsilon * upper_bound))
        
        # Add contribution for x values inside the current interval
        idx = (x_grid > lower_bound) & (x_grid <= upper_bound)
        integral_N_term[idx] += (count / epsilon) * (np.exp(-epsilon * lower_bound) - np.exp(-epsilon * x_grid[idx]))

    integral_total = integral_C_term + epsilon**2 * integral_N_term
    
    # Combine terms to get the final solution y(x)
    y = np.exp(epsilon * x_grid) * (1 + integral_total)
    
    # Enforce boundary conditions to correct for any minor numerical drift
    y[0] = 1.0
    y[-1] = 0.0
    return y

if __name__ == '__main__':
    # Run the main function to get the scaling exponent
    scaling_exponent = estimate_fluctuation_scaling()
