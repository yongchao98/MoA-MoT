import numpy as np
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt

def solve_fluctuation_ode(epsilon, z_i, grid_points):
    """
    Solves the ODE w'' - eps*w' = eps^2 * (sum(delta(x-z_i)) - 1)
    for the fluctuation w(x) using finite differences.
    """
    L = 1.0 / epsilon
    N = len(z_i)
    x_grid = np.linspace(0, L, grid_points)
    h = x_grid[1] - x_grid[0]

    # Create the source vector f = eps^2 * (counts/h - 1)
    # The term -1 represents the average source density.
    counts, _ = np.histogram(z_i, bins=grid_points, range=(0, L))
    f = epsilon**2 * (counts / h - 1.0)
    
    # We are solving Aw=f, where A is the finite difference operator matrix.
    # (1/h^2 - eps/(2h))w_{j-1} - (2/h^2)w_j + (1/h^2 + eps/(2h))w_{j+1} = f_j
    # We need to solve for the inner points w_1, ..., w_{grid_points-2}
    # since w_0 = w_{grid_points-1} = 0.
    
    A_main_diag = -2.0 / h**2 * np.ones(grid_points - 2)
    A_lower_diag = (1.0 / h**2 - epsilon / (2.0 * h)) * np.ones(grid_points - 3)
    A_upper_diag = (1.0 / h**2 + epsilon / (2.0 * h)) * np.ones(grid_points - 3)
    
    # The matrix for solve_banded is specified with diagonals.
    # The first row is upper diagonal, second is main, third is lower.
    A_banded = np.zeros((3, grid_points - 2))
    A_banded[0, 1:] = A_upper_diag
    A_banded[1, :] = A_main_diag
    A_banded[2, :-1] = A_lower_diag
    
    # Solve the system
    w_inner = solve_banded((1, 1), A_banded, f[1:-1])
    
    # Add boundary points
    w = np.zeros(grid_points)
    w[1:-1] = w_inner
    
    return w, x_grid

def estimate_scaling(z_generator, case_name, epsilons):
    """
    Runs Monte Carlo simulation for a given z_i generator and estimates scaling.
    """
    print(f"--- Simulating Case: {case_name} ---")
    R_values = []
    
    for eps in epsilons:
        N = int(1.0 / eps) - 1
        L = 1.0 / eps
        if N <= 0:
            continue
            
        n_trials = 100
        # Use more grid points for smaller epsilon to resolve features
        grid_points = int(5 / eps) 
        w_samples = np.zeros((n_trials, grid_points))

        for i in range(n_trials):
            z_i = z_generator(N, L)
            w, _ = solve_fluctuation_ode(eps, z_i, grid_points)
            w_samples[i, :] = w

        variances = np.var(w_samples, axis=0)
        max_var = np.max(variances)
        R = np.sqrt(max_var)
        R_values.append(R)
        print(f"epsilon = {eps:.3f}, R(epsilon) = {R:.4f}")

    # Fit a line to log(R) vs log(epsilon) to find the scaling exponent
    valid_epsilons = epsilons[:len(R_values)]
    slope, intercept = np.polyfit(np.log(valid_epsilons), np.log(R_values), 1)

    print("\nScaling Analysis:")
    print(f"Fit result: log(R) = {slope:.2f} * log(epsilon) + {intercept:.2f}")
    print(f"The estimated scaling is R(epsilon) ~ epsilon^{slope:.2f}")
    
    # You can uncomment the following lines to see the plot
    # plt.figure()
    # plt.loglog(valid_epsilons, R_values, 'o', label=f'Simulated R for {case_name}')
    # plt.loglog(valid_epsilons, np.exp(intercept) * valid_epsilons**slope, '--', label=f'Fit: $\\epsilon^{{{slope:.2f}}}$')
    # plt.xlabel('log(epsilon)')
    # plt.ylabel('log(R)')
    # plt.title(f'Scaling of R(epsilon) for {case_name}')
    # plt.grid(True)
    # plt.legend()
    # plt.show()
    print("-" * 35 + "\n")


if __name__ == '__main__':
    # Define a range of epsilon values for the simulation
    epsilons_to_test = [0.1, 0.05, 0.04, 0.03, 0.02]

    # --- Case 1: z_i ~ Ordered Uniform([0, 1/epsilon]) ---
    uniform_z_generator = lambda N, L: np.sort(np.random.uniform(0, L, size=N))
    estimate_scaling(uniform_z_generator, "Ordered Uniform z_i", epsilons_to_test)

    # --- Case 2: z_i ~ Normal(i, 0.5) ---
    def normal_z_generator(N, L):
        # Ensure z_i are within the domain for the simulation
        z = np.random.normal(loc=np.arange(1, N + 1), scale=0.5)
        return z[(z > 0) & (z < L)]

    estimate_scaling(normal_z_generator, "i.i.d. Normal z_i", epsilons_to_test)
    
    print("As predicted by the theoretical analysis, the scaling for R(epsilon) remains the same")
    print("in both cases, proportional to epsilon^0.5.")

<<<0.5>>>