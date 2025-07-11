import numpy as np

def calculate_variance_of_y():
    """
    Calculates the variance of Y through Monte Carlo simulation.

    Let X_1, X_2, X_3, and X_4 be i.i.d. random variables U(0, 1).
    Y is the second closest value to X_1 among X_2, X_3, and X_4.
    Var(Y) = E[Y^2] - (E[Y])^2.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Number of simulation trials
    num_simulations = 5_000_000

    # Generate all random numbers at once for efficiency
    # samples is a (num_simulations x 4) array
    samples = np.random.rand(num_simulations, 4)

    # Extract x1 and the other three points (x2, x3, x4)
    x1 = samples[:, 0:1]  # Shape: (N, 1) for broadcasting
    others = samples[:, 1:] # Shape: (N, 3)

    # Calculate the absolute distances from x1 to each of the other points
    # distances will have shape (N, 3)
    distances = np.abs(others - x1)

    # For each row, get the indices that would sort the distances
    # E.g., for a row of distances [0.3, 0.1, 0.4], argsort gives [1, 0, 2]
    sorted_dist_indices = np.argsort(distances, axis=1)

    # The index of the second closest point is the middle column of the sorted indices
    second_closest_idx = sorted_dist_indices[:, 1]

    # Use the indices to select the corresponding value of Y from the 'others' array.
    # This uses advanced NumPy indexing.
    # np.arange(num_simulations) provides the row index for each simulation.
    # second_closest_idx provides the column index (0, 1, or 2) for the 'others' array.
    y_values = others[np.arange(num_simulations), second_closest_idx]

    # Calculate the moments and variance from the simulated y_values
    mean_y = np.mean(y_values)
    mean_y_squared = np.mean(y_values**2)
    
    # We can calculate variance from the moments or using np.var
    # var_y = mean_y_squared - mean_y**2
    var_y = np.var(y_values) # ddof=0 for population variance, which is what we estimate

    # Output the results based on the variance equation
    print("Based on the simulation with {} trials:".format(num_simulations))
    print("Estimated E[Y] = {:.6f}".format(mean_y))
    print("Estimated E[Y^2] = {:.6f}".format(mean_y_squared))
    print("\nThe variance is calculated as Var(Y) = E[Y^2] - (E[Y])^2")
    print("Var(Y) = {:.6f} - ({:.6f})^2 = {:.6f}".format(mean_y_squared, mean_y, var_y))
    print("\nThe estimated variance of Y is {:.6f}".format(var_y))
    # The analytical result is 7/60.
    # print("The exact analytical result is 7/60 ≈ {:.6f}".format(7/60))

calculate_variance_of_y()