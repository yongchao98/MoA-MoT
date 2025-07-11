import numpy as np

def estimate_variance():
    """
    This function estimates the variance of Y using a Monte Carlo simulation.
    
    Y is defined as the second closest value to X1 among X2, X3, and X4,
    where all X_i are i.i.d. random variables from U(0, 1).
    """
    # Step 1: Set the number of simulations for accuracy.
    num_simulations = 2000000

    # Step 2: Generate all random numbers in a single batch for efficiency.
    # Each row represents a set of (X1, X2, X3, X4).
    samples = np.random.rand(num_simulations, 4)

    # Separate X1 from the other variables.
    # We use slicing to keep the dimensions correct for broadcasting.
    x1 = samples[:, 0:1]
    x_others = samples[:, 1:]

    # Step 3: Calculate the distances of X2, X3, X4 from X1 for each simulation.
    distances = np.abs(x_others - x1)

    # Step 4: For each simulation, find the indices that would sort the distances.
    # The index of the second smallest distance will be at the 2nd position (index 1).
    sorted_indices = np.argsort(distances, axis=1)
    second_closest_idx = sorted_indices[:, 1]

    # Step 5: Use the indices to retrieve the value of Y for each simulation.
    # np.arange creates row indices [0, 1, ..., N-1] to pair with our column indices.
    y_values = x_others[np.arange(num_simulations), second_closest_idx]

    # Step 6: Calculate the variance of the collected Y values.
    # We calculate the mean of Y and the mean of Y^2 to show the full variance equation.
    mean_y = np.mean(y_values)
    mean_y_squared = np.mean(y_values**2)
    variance_y = np.var(y_values) # np.var is numerically stable

    # Step 7: Output the results as requested.
    print(f"The variance is calculated using the formula: Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Based on {num_simulations} simulations:")
    print(f"  - Estimated E[Y] = {mean_y:.6f}")
    print(f"  - Estimated E[Y^2] = {mean_y_squared:.6f}")
    
    # Print the final equation with the calculated numbers.
    print("\nFinal Equation:")
    print(f"{variance_y:.6f} = {mean_y_squared:.6f} - ({mean_y:.6f})^2")

    # The final answer. The theoretical value is 1/30.
    print("\nFinal estimated variance:")
    print(variance_y)


if __name__ == '__main__':
    estimate_variance()