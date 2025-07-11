import numpy as np

def count_kk_modes(threshold, num_points=1000):
    """
    Calculates the number of Kaluza-Klein spin 2 modes with mass-squared
    below a given threshold.

    Args:
        threshold (float): The mass-squared threshold.
        num_points (int): The number of grid points for the numerical solution.
    
    Returns:
        int: The number of eigenvalues below the threshold.
    """
    # 1. Set up the numerical grid for the extra dimension x in [0, 2*pi]
    # We use endpoint=False because the functions are periodic.
    # dx is the step size.
    x_grid, dx = np.linspace(0, 2 * np.pi, num_points, endpoint=False, retstep=True)

    # We use a 'staggered' grid for the function p(x), evaluating it at
    # points halfway between the main grid points (x_{j+1/2}).
    x_half_step_grid = x_grid + dx / 2.0

    # 2. Define the functions from the problem statement.
    # The warp factor is A(x).
    def warp_factor_A(x):
        return np.sin(x) + 4 * np.cos(x)

    # The function p(x) in the differential operator.
    def p_function(x):
        return np.exp(2 * warp_factor_A(x))

    # Evaluate p(x) on the half-step grid.
    p_half = p_function(x_half_step_grid)
    
    # 3. Construct the finite difference matrix M for the operator -d/dx(p(x)d/dx).
    # This matrix represents the operator on the discretized grid.
    M = np.zeros((num_points, num_points))
    dx_squared = dx * dx

    # The periodic boundary conditions are handled using the modulo operator (%)
    # on the indices.
    for i in range(num_points):
        # Diagonal elements M_ii
        # This corresponds to the part of the operator acting on psi_i.
        p_plus_half = p_half[i]
        p_minus_half = p_half[i - 1] # Relies on Python's negative indexing for i=0
        M[i, i] = (p_plus_half + p_minus_half) / dx_squared
        
        # Off-diagonal elements M_{i, i+1} and M_{i, i-1}
        # These correspond to the parts of the operator acting on psi_{i+1} and psi_{i-1}.
        M[i, (i + 1) % num_points] = -p_plus_half / dx_squared
        M[i, (i - 1 + num_points) % num_points] = -p_minus_half / dx_squared

    # 4. Calculate the eigenvalues of the matrix M.
    # The eigenvalues of M are the approximate values of m^2.
    # We use eigvalsh because the matrix is real and symmetric (Hermitian), which
    # is faster and more numerically stable.
    eigenvalues = np.linalg.eigvalsh(M)

    # 5. Count the number of eigenvalues below the threshold.
    count = np.sum(eigenvalues < threshold)
    
    return count

# Set the threshold for the mass-squared
mass_squared_threshold = 14

# Calculate the number of modes
num_modes = count_kk_modes(mass_squared_threshold)

# Print the final result
print(f"The warp factor is given by A(x) = sin(x) + 4*cos(x).")
print(f"We need to find the number of eigenvalues (m^2) below the threshold value of {mass_squared_threshold}.")
print(f"Number of eigenvalues below {mass_squared_threshold}: {num_modes}")
