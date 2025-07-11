import numpy as np
from scipy.stats import linregress
import warnings

# Suppress RankWarning from polyfit which can happen with few data points
warnings.simplefilter('ignore', np.RankWarning)

def estimate_R_scaling():
    """
    Performs a Monte Carlo simulation to estimate the scaling of R(epsilon)
    for the uniform distribution case.

    R(epsilon) = (max_x |Var[y(x) - y(0)]|)^1/2

    The function simulates the solution for several epsilon values, calculates R,
    and then fits a power law R(epsilon) = C * epsilon^p to the results.
    """
    
    # We choose epsilon values on a log scale for better fitting of the power law.
    # We avoid very small epsilons which would make N too large and the simulation slow.
    epsilon_values = np.array([0.2, 0.1, 0.05, 0.04])
    mc_runs = 250  # Number of Monte Carlo runs for each epsilon

    R_values = []
    actual_epsilons = []

    print("Thinking Process: Running Monte Carlo simulations...")
    for epsilon in epsilon_values:
        L = 1.0 / epsilon
        N = int(np.floor(L - 1))
        
        if N <= 1:
            print(f"  eps = {epsilon:.3f} is too large (N={N}), skipping.")
            continue
        
        print(f"  Simulating for eps = {epsilon:.3f} (L = {L:.1f}, N = {N}) with {mc_runs} runs...")
        
        # Grid of x-points to evaluate the solution and its variance
        x_grid = np.linspace(0, L, 250)
        y_samples = []

        for _ in range(mc_runs):
            # For each run, generate new random source locations
            z = np.sort(np.random.uniform(0, L, N))
            
            # --- Set up and solve the linear system for coefficients ---
            # In each of the N+1 intervals [z_j, z_{j+1}], y(x) = A_j + B_j * exp(eps*x).
            # This requires solving a 2*(N+1) x 2*(N+1) linear system.
            num_coeffs = 2 * (N + 1)
            M = np.zeros((num_coeffs, num_coeffs))
            b = np.zeros(num_coeffs)

            # Boundary Condition at x=0: y_0(0) = 1
            M[0, 0] = 1
            M[0, 1] = 1
            b[0] = 1

            # z_i are the N internal boundaries
            for i in range(N):
                z_i = z[i]
                exp_term = np.exp(epsilon * z_i)
                # Continuity y_i(z_i) = y_{i+1}(z_i)
                row1 = 2 * i + 1
                M[row1, 2*i:2*i+2] = [1, exp_term]
                M[row1, 2*(i+1):2*(i+1)+2] = [-1, -exp_term]
                # Derivative jump y'_{i+1}(z_i) - y'_i(z_i) = eps^2
                row2 = 2 * i + 2
                M[row2, 2*i + 1] = -epsilon * exp_term
                M[row2, 2*(i+1) + 1] = epsilon * exp_term
                b[row2] = epsilon**2

            # Boundary Condition at x=L: y_N(L) = 0
            M[-1, -2] = 1
            M[-1, -1] = np.exp(epsilon * L)

            try:
                coeffs = np.linalg.solve(M, b)
            except np.linalg.LinAlgError:
                continue # Skip run if matrix is singular

            # --- Evaluate the solution on the grid ---
            y_vals = np.zeros_like(x_grid)
            z_extended = np.concatenate(([0], z, [L]))
            
            # This loop is faster than calling searchsorted for each point individually
            for j in range(N + 1):
                # Find all x that fall in interval j
                mask = (x_grid >= z_extended[j]) & (x_grid <= z_extended[j+1])
                A_j, B_j = coeffs[2*j], coeffs[2*j+1]
                y_vals[mask] = A_j + B_j * np.exp(epsilon * x_grid[mask])
            
            y_samples.append(y_vals)

        if not y_samples:
            print(f"    Warning: all runs failed for epsilon = {epsilon}")
            continue
            
        y_samples = np.array(y_samples)
        
        # y(0) is a constant (1), so Var[y(x) - y(0)] = Var[y(x)].
        variances = np.var(y_samples, axis=0)
        max_variance = np.max(variances)
        
        R = np.sqrt(max_variance)
        R_values.append(R)
        actual_epsilons.append(epsilon)
        print(f"    -> Result: R({epsilon:.3f}) = {R:.4f}")

    # --- Fit the results to a power law R = C * epsilon^p ---
    # This is equivalent to a linear fit on log-log data: log(R) = p*log(eps) + log(C)
    log_R = np.log(R_values)
    log_eps = np.log(actual_epsilons)
    
    # Use polyfit to get slope (p) and intercept (log(C))
    p, log_C = np.polyfit(log_eps, log_R, 1)
    C = np.exp(log_C)
    
    print("\n--- Final Scaling Result ---")
    print(f"Numerical simulation suggests the scaling law is R(epsilon) = C * epsilon^p")
    print("The estimated parameters are:")
    print(f"  C = {C:.4f}")
    print(f"  p = {p:.4f}")

    print("\nThe estimated maximum magnitude of fluctuations as a function of epsilon is:")
    final_equation = f"R(epsilon) = {C:.4f} * epsilon^{p:.4f}"
    print(final_equation)


if __name__ == '__main__':
    estimate_R_scaling()
<<<0.5>>>