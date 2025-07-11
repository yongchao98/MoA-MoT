import numpy as np

def solve_ode_coefficients(epsilon, z_coords):
    """
    Calculates the coefficients A_k, B_k for the piecewise solution
    y_k(x) = A_k + B_k * exp(epsilon * x).
    
    The solution is constructed by enforcing continuity of y and a specific
    jump in y' at each source location z_i. The coefficients are uniquely
    determined by the boundary conditions.
    """
    N = len(z_coords)
    e = np.e

    # Calculate sum(exp(eps*z_i)) for i=1 to N
    sum_exp_z = np.sum(np.exp(epsilon * z_coords))

    # Solve for the initial coefficient B_0 using the boundary condition at x=1/epsilon.
    # The analytical derivation gives: B0(e-1) = -1 + eps * sum(exp(eps*z_i)) - N*eps*e
    B0 = (-1.0 + epsilon * sum_exp_z - N * epsilon * e) / (e - 1.0)
    A0 = 1.0 - B0

    # The coefficients for subsequent intervals are found by recurrence relations.
    A_coeffs = [A0]
    B_coeffs = [B0]
    current_A = A0
    current_B = B0
    for i in range(N):
        # Jump conditions from integrating the ODE across a delta function: 
        # y is continuous, y' jumps by epsilon^2 at each z_i.
        # This leads to: B_k = B_{k-1} + epsilon
        #                A_k = A_{k-1} - epsilon * exp(epsilon * z_i)
        current_B = current_B + epsilon
        current_A = current_A - epsilon * np.exp(epsilon * z_coords[i])
        A_coeffs.append(current_A)
        B_coeffs.append(current_B)
        
    return A_coeffs, B_coeffs

def get_solution_function(epsilon, z_coords, A_coeffs, B_coeffs):
    """
    Returns a function that computes y(x) for a given vector of x values
    using the pre-calculated coefficients.
    """
    def y_solution(x_vals):
        # Ensure x_vals is a numpy array for vectorized operations
        if not isinstance(x_vals, np.ndarray):
            x_vals = np.array([x_vals])
        
        y_vals = np.zeros_like(x_vals, dtype=float)
        
        # `searchsorted` finds the interval index for each x value.
        # If x is in (z_{i-1}, z_i), the index is i.
        # The z_coords must be sorted for this to work.
        indices = np.searchsorted(z_coords, x_vals)
        
        # Calculate y(x) for each interval using the corresponding coefficients
        for i, (A, B) in enumerate(zip(A_coeffs, B_coeffs)):
            mask = (indices == i)
            if np.any(mask):
                y_vals[mask] = A + B * np.exp(epsilon * x_vals[mask])
        return y_vals

    return y_solution

def estimate_fluctuation_scaling():
    """
    Main function to perform the Monte Carlo simulation and estimate the scaling of R.
    """
    # Epsilon values to test. They should be small and span an order of magnitude.
    epsilons = [0.08, 0.04, 0.02, 0.01]
    # Number of Monte Carlo trials to average over for each epsilon
    num_mc_trials = 100 
    # Number of spatial points to evaluate the solution on
    num_x_points = 1000   
    
    results_R = []

    print("Starting numerical estimation of the scaling of R with epsilon...")
    
    for eps in epsilons:
        N = int(1.0/eps) - 1
        L = 1.0/eps
        x_grid = np.linspace(0, L, num_x_points)
        y0_grid = 1.0 - eps * x_grid
        
        # Accumulators for calculating variance: Var(X) = E[X^2] - (E[X])^2
        fluctuations_sq_sum = np.zeros(num_x_points)
        fluctuations_sum = np.zeros(num_x_points)

        print(f"\nRunning for epsilon = {eps} (N={N}, L={L:.1f})")
        for trial in range(num_mc_trials):
            # Generate N ordered uniform random numbers for the source locations
            z_i = np.sort(np.random.uniform(0, L, N))
            
            # Solve the ODE for this specific realization of z_i
            A_coeffs, B_coeffs = solve_ode_coefficients(eps, z_i)
            y_func = get_solution_function(eps, z_i, A_coeffs, B_coeffs)
            y_vals = y_func(x_grid)
            
            # Calculate fluctuation y(x) - y0(x)
            fluc = y_vals - y0_grid
            
            # Accumulate sum and sum of squares
            fluctuations_sum += fluc
            fluctuations_sq_sum += fluc**2

        # Calculate variance across the Monte Carlo trials
        mean_fluc = fluctuations_sum / num_mc_trials
        mean_sq_fluc = fluctuations_sq_sum / num_mc_trials
        variance = np.abs(mean_sq_fluc - mean_fluc**2)
        
        # Estimate R = sqrt(max(Var))
        max_var = np.max(variance)
        R = np.sqrt(max_var)
        results_R.append(R)
        print(f"  --> Estimated R = {R:.5f}")

    # --- Analysis of Results ---
    # Fit a line to the log-log plot to find the scaling exponent
    log_eps = np.log(epsilons)
    log_R = np.log(results_R)
    
    # Fit line: log(R) = p*log(eps) + log(C)
    p, log_C = np.polyfit(log_eps, log_R, 1)
    C = np.exp(log_C)
    
    print("\n-------------------------------------------------------------")
    print("Log-log regression analysis to find the scaling R = C * epsilon^p")
    print("-------------------------------------------------------------")
    print(f"The final equation for the fluctuation magnitude is R(epsilon) = C * epsilon^p")
    print(f"The estimated value for the coefficient C is: {C:.4f}")
    print(f"The estimated value for the scaling exponent p is: {p:.4f}")
    print("\nTheoretical analysis suggests the scaling exponent should be p = 0.5.")
    
    print("\n--- Answer to the second part of the question ---")
    print("Q: Do you expect the scaling for R(epsilon) to remain the same if z_i ~ Normal(i, 0.5)?")
    print("A: Yes, I expect the scaling to remain the same. The magnitude of the fluctuations is primarily determined by the total number of random sources, N, which scales as epsilon^-1 in both cases. The variance of the solution scales proportionally to N. This leads to Var[y] ~ epsilon, and therefore R ~ epsilon^0.5. The specific distribution of the sources (uniform vs. a train of Gaussians) would change the pre-factor C, but not the fundamental scaling exponent p.")
    print("-------------------------------------------------------------\n")

if __name__ == '__main__':
    estimate_fluctuation_scaling()