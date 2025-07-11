import numpy as np

def solve_one_realization(x_grid, z, epsilon):
    """
    Calculates the exact solution y(x) for a given realization of source locations z.
    
    The solution is y(x) = C1 + C2*exp(eps*x) + eps * Sum_{i} [exp(eps*(x-z_i))-1]*H(x-z_i),
    where H is the Heaviside step function and C1, C2 are constants determined
    by the boundary conditions y(0)=1, y(1/eps)=0.
    """
    L = 1.0 / epsilon
    e = np.e # Euler's number
    
    # Calculate the sum term needed for the integration constants
    # The sum is over all sources z_i
    sum_exp_z = np.sum(np.exp(-epsilon * z))
    
    # Calculate constants C1 and C2 from boundary conditions
    C2 = -epsilon * (1.0 + e * sum_exp_z) / (e - 1.0)
    C1 = 1.0 - C2
    
    y = np.zeros_like(x_grid)
    
    # Calculate y(x) for each point in the x_grid
    for i, x in enumerate(x_grid):
        # Homogeneous part of the solution
        y_h = C1 + C2 * np.exp(epsilon * x)
        
        # Particular solution part (sum over sources z_i < x)
        active_z = z[z < x]
        if len(active_z) > 0:
            sum_term = np.sum(np.exp(epsilon * (x - active_z)) - 1.0)
        else:
            sum_term = 0.0
            
        y_p = epsilon * sum_term
        
        y[i] = y_h + y_p
        
    return y

def estimate_scaling_exponent():
    """
    Main function to run the simulation for different epsilon values,
    calculate R, and fit a scaling law R = C * epsilon^p.
    """
    # Use a fixed random seed for reproducibility
    np.random.seed(42)

    # A set of epsilon values to test
    epsilons = np.array([0.1, 0.08, 0.05, 0.04, 0.025])
    
    # Simulation parameters
    num_realizations = 2000
    num_x_points = 1000
    
    results_R = []

    print("Starting simulation...")
    for eps in epsilons:
        L = 1.0 / eps
        # Ensure N is an integer
        N = int(round(1.0 / eps) - 1.0)
        if N <= 0:
            continue

        print(f"Running for epsilon = {eps:.3f} (N={N}, L={L:.1f})")
        
        x_grid = np.linspace(0, L, num_x_points)
        y_realizations = np.zeros((num_realizations, num_x_points))
        y0_on_grid = 1.0 - eps * x_grid
        
        for i in range(num_realizations):
            # Generate N i.i.d. random variables from Uniform(0, L) and sort them
            z_iid = np.random.uniform(0, L, N)
            z = np.sort(z_iid)
            
            y = solve_one_realization(x_grid, z, eps)
            y_realizations[i, :] = y
        
        # Calculate Var[y(x) - y0(x)] which is equal to Var[y(x)]
        # The variance is calculated along axis 0 (across realizations)
        var_y = np.var(y_realizations, axis=0)
        
        # R is the square root of the maximum variance
        R = np.sqrt(np.max(var_y))
        results_R.append(R)
        print(f"-> Calculated R = {R:.5f}")

    # Fit a line to the log-log data to find the scaling exponent
    log_eps = np.log(epsilons)
    log_R = np.log(results_R)
    
    # polyfit returns [slope, intercept] for a 1-degree polynomial
    coeffs = np.polyfit(log_eps, log_R, 1)
    
    scaling_exponent = coeffs[0]
    prefactor_C = np.exp(coeffs[1])

    print("\n--- Simulation Complete ---")
    print("The scaling of the fluctuation magnitude R with epsilon is estimated.")
    print("The relationship is modeled as: R = C * epsilon^p")
    print("\nFitted Parameters:")
    print(f"Prefactor C = {prefactor_C:.4f}")
    print(f"Scaling Exponent p = {scaling_exponent:.4f}")
    print("\nFinal Equation:")
    print(f"R(epsilon) ≈ {prefactor_C:.4f} * epsilon^{scaling_exponent:.4f}")
    
# Execute the estimation
estimate_scaling_exponent()