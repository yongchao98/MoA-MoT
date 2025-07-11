import numpy as np
from scipy.linalg import solve_banded
from scipy.stats import linregress
import matplotlib.pyplot as plt

def solve_bvp(epsilon, z_points):
    """
    Solves the BVP y'' - epsilon*y' = f(x) using finite differences.
    f(x) is a sum of delta functions.
    """
    L = 1.0 / epsilon
    N = int(L) - 1
    
    # Finite difference grid
    grid_size = 2000 # Number of grid points
    x_grid = np.linspace(0, L, grid_size)
    h = x_grid[1] - x_grid[0]
    
    # Construct the RHS vector f(x)
    rhs_f = np.zeros(grid_size)
    # Discretize the delta functions
    z_indices = np.floor(np.array(z_points) / h).astype(int)
    # Ensure indices are within bounds
    z_indices = z_indices[z_indices < grid_size]
    counts = np.bincount(z_indices, minlength=grid_size)
    # The magnitude of f(x) is epsilon^2 * sum(delta)
    # The density is counts/h
    rhs_f = epsilon**2 * counts / h

    # Construct the tridiagonal matrix for the linear system Ay=b
    # (1/h^2 - e/2h)y_{j-1} + (-2/h^2)y_j + (1/h^2 + e/2h)y_{j+1} = f_j
    # We solve for the internal points y[1:-1]
    num_internal_points = grid_size - 2
    
    # Create banded matrix diagonals
    lower_diag = np.full(num_internal_points, 1/h**2 - epsilon/(2*h))
    main_diag = np.full(num_internal_points, -2/h**2)
    upper_diag = np.full(num_internal_points, 1/h**2 + epsilon/(2*h))
    
    ab = np.vstack([np.pad(upper_diag, (1,0), 'constant'), 
                    main_diag, 
                    np.pad(lower_diag, (0,1), 'constant')])

    # Adjust RHS for boundary conditions y(0)=1, y(L)=0
    rhs_b = rhs_f[1:-1]
    # y(0) = 1
    rhs_b[0] -= lower_diag[0] * 1.0 
    # y(L) = 0
    # rhs_b[-1] -= upper_diag[-1] * 0.0 (no change needed)

    # Solve the banded system
    y_internal = solve_banded((1, 1), ab, rhs_b)
    
    # Combine with boundary points to get full solution
    y_solution = np.hstack([1.0, y_internal, 0.0])
    
    return x_grid, y_solution

def estimate_R_for_epsilon(epsilon, num_simulations=100):
    """
    Runs Monte Carlo simulations for a given epsilon to estimate R.
    """
    L = 1.0 / epsilon
    N = int(L) - 1
    grid_size = 2000
    
    y_samples = np.zeros((num_simulations, grid_size))
    
    for i in range(num_simulations):
        # Generate N ordered uniform random variables
        z_i = np.sort(np.random.uniform(0, L, N))
        
        x, y = solve_bvp(epsilon, z_i)
        y_samples[i, :] = y
        
    # Var[y(x) - y(0)] = Var[y(x)] since y(0)=1 is constant.
    var_y = np.var(y_samples, axis=0)
    max_var = np.max(var_y)
    R = np.sqrt(max_var)
    return R

def main():
    """
    Main function to run the simulation and find the scaling law.
    """
    # Epsilon values to test (smaller values require more computation)
    epsilons = [0.1, 0.08, 0.05, 0.04, 0.03]
    R_values = []
    
    print("Running simulations for different epsilon values...")
    for eps in epsilons:
        r_val = estimate_R_for_epsilon(eps, num_simulations=100)
        R_values.append(r_val)
        print(f"  epsilon = {eps:.3f}, R(epsilon) = {r_val:.5f}")

    # Fit a power law R = C * epsilon^a using log-log regression
    log_eps = np.log(epsilons)
    log_R = np.log(R_values)
    slope, intercept, _, _, _ = linregress(log_eps, log_R)
    
    a = slope
    C = np.exp(intercept)

    print("\n--- Results ---")
    print("The relationship between R and epsilon is modeled as R = C * epsilon^a")
    print(f"Theoretical prediction: R = (1/sqrt(192)) * epsilon^0.5 ≈ 0.072 * epsilon^0.5")

    print("\nFitted power law from numerical simulation:")
    # The final code should output each number in the final equation.
    print(f"Final Equation: R(epsilon) = {C:.4f} * epsilon^{a:.4f}")

    # Plotting for visualization
    plt.figure(figsize=(8, 6))
    plt.scatter(log_eps, log_R, label='Simulation Results', color='blue')
    plt.plot(log_eps, intercept + slope * log_eps, 'r-', label=f'Fit: R={C:.3f}*eps^{a:.3f}')
    plt.plot(log_eps, np.log(1/np.sqrt(192)) + 0.5 * log_eps, 'g--', label='Theoretical: R ~ eps^0.5')
    plt.xlabel('log(epsilon)')
    plt.ylabel('log(R)')
    plt.title('Log-Log plot of R vs. epsilon')
    plt.legend()
    plt.grid(True)
    # plt.show() # Uncomment to display the plot

if __name__ == '__main__':
    main()
