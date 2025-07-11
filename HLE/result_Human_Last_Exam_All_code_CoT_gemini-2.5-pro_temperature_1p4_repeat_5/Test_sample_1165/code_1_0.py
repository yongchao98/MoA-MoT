import numpy as np
from scipy.stats import linregress

def greens_function_solution(x, L, z):
    """
    Calculates the solution Y(x) = sum_i G_0(x, z_i) for a given x.
    G_0 is the Green's function for Y''=f with Y(0)=Y(L)=0.
    """
    z_less = z[z < x]
    z_more_equal = z[z >= x]
    
    # G_0(x,s) = s(x-L)/L for s < x
    sum1 = np.sum(z_less * (x - L) / L)
    
    # G_0(x,s) = x(s-L)/L for s >= x
    sum2 = np.sum(x * (z_more_equal - L) / L)
    
    return sum1 + sum2

def estimate_R_for_epsilon(epsilon, num_mc_runs=1000, num_x_points=200):
    """
    Estimates the fluctuation magnitude R for a given epsilon using Monte Carlo simulation.
    """
    if epsilon >= 1 or epsilon <= 0:
        raise ValueError("Epsilon must be between 0 and 1.")
        
    L = 1.0 / epsilon
    # As per problem, N = epsilon^-1 - 1, which is L-1
    N = int(np.floor(L - 1))
    
    if N <= 0:
        # For large epsilon, N can be 0 or negative.
        # The model is not well-defined here, so we return NaN.
        return np.nan

    x_grid = np.linspace(0, L, num_x_points)
    
    # Each row will be a single realization of Y(x) over the x_grid
    y2_stoch_samples = np.zeros((num_mc_runs, num_x_points))

    for i in range(num_mc_runs):
        # Generate N ordered values from Uniform([0, L])
        # This is done by generating N i.i.d samples and sorting them.
        z_i = np.sort(np.random.uniform(0, L, N))
        
        for j, x_val in enumerate(x_grid):
            # The stochastic part of y2 solves Y'' = sum(delta(x-z_i))
            y2_stoch_samples[i, j] = greens_function_solution(x_val, L, z_i)
            
    # Calculate variance at each point x across all Monte Carlo runs
    var_y2 = np.var(y2_stoch_samples, axis=0)
    
    # Find the maximum variance
    max_var_y2 = np.max(var_y2)
    
    # R = epsilon^2 * sqrt(max(Var[y2]))
    R = epsilon**2 * np.sqrt(max_var_y2)
    
    return R

def main():
    """
    Runs the simulation for a range of epsilon values and determines the scaling law.
    """
    # A range of epsilon values to test, on a log scale
    epsilons = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
    
    print("Running simulations to estimate R(epsilon)...")
    
    estimated_Rs = []
    for eps in epsilons:
        R = estimate_R_for_epsilon(eps, num_mc_runs=500, num_x_points=200)
        estimated_Rs.append(R)
        print(f"  For epsilon = {eps:.5f}, estimated R = {R:.5f}")
        
    estimated_Rs = np.array(estimated_Rs)
    
    # Filter out any NaN results that might occur for larger epsilon
    valid_indices = ~np.isnan(estimated_Rs)
    if np.sum(valid_indices) < 2:
        print("\nNot enough data points to determine scaling law.")
        return

    valid_epsilons = epsilons[valid_indices]
    valid_Rs = estimated_Rs[valid_indices]

    # Perform a linear regression on the log-log data to find the scaling exponent
    log_eps = np.log(valid_epsilons)
    log_R = np.log(valid_Rs)
    
    slope, intercept, _, _, _ = linregress(log_eps, log_R)
    
    # The scaling is R = C * epsilon^p
    # log(R) = log(C) + p * log(epsilon)
    # p is the slope, C = exp(intercept)
    p = slope
    C = np.exp(intercept)

    print("\n--- Results ---")
    print("The relationship is estimated to be of the form: R(epsilon) = C * epsilon^p")
    print("\nFinal Equation:")
    print(f"R(epsilon) = {C:.4f} * epsilon^{p:.4f}")
    
    # Analytical prediction
    # R = epsilon^(1/2) / sqrt(192)
    # C_analytical = 1 / sqrt(192) = 1 / (8*sqrt(3)) approx 0.072
    # p_analytical = 0.5
    print("\nAnalytical prediction for comparison:")
    C_analytical = 1 / np.sqrt(192)
    p_analytical = 0.5
    print(f"R(epsilon) = {C_analytical:.4f} * epsilon^{p_analytical:.4f}")


    print("\n--- Answer to the second question ---")
    print("If z_i were drawn from Normal(i, 0.5), the scaling for R(epsilon) would NOT remain the same.")
    print("The increased regularity of the source locations z_i would reduce the overall variance.")
    print("Analytical calculations show that the scaling would change from epsilon^0.5 to epsilon^1.5, meaning the fluctuations would be significantly smaller.")

if __name__ == '__main__':
    main()
