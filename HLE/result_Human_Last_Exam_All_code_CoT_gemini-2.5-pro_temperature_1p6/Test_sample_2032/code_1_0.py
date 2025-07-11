import numpy as np

def find_variance_of_y():
    """
    This function calculates the variance of Y using a Monte Carlo simulation.

    Let X_1, X_2, X_3, and X_4 be i.i.d. random variables uniformly
    distributed on the interval [0, 1]. Y is the second closest value to
    X_1 among X_2, X_3, and X_4.
    """
    # Number of simulations. A large number provides a more accurate estimate.
    num_simulations = 10_000_000

    print(f"Running Monte Carlo simulation with {num_simulations:,} samples...")

    # Step 1: Generate N samples for X1, X2, X3, X4 from U(0,1)
    # This creates an array of shape (num_simulations, 4)
    samples = np.random.rand(num_simulations, 4)

    # Step 2: Separate X1 from the other variables (X2, X3, X4)
    x1 = samples[:, 0]
    x_others = samples[:, 1:]

    # Step 3: Calculate the absolute distances from X1 for each sample
    # x1 is shape (N,), x_others is shape (N, 3). We reshape x1 to (N, 1)
    # to allow numpy's broadcasting to subtract it from each of the 3 columns.
    distances = np.abs(x_others - x1[:, np.newaxis])

    # Step 4: For each sample, find the index of the second-smallest distance.
    # np.argsort returns the indices that would sort the array.
    # axis=1 sorts along each row (i.e., for each of the N samples).
    # The indices will be 0, 1, or 2, corresponding to X2, X3, or X4.
    sorted_indices = np.argsort(distances, axis=1)
    
    # The index of the second-smallest distance is in the second column (index 1).
    second_closest_idx = sorted_indices[:, 1]

    # Step 5: Use these indices to select the value of Y for each sample.
    # We create an array [0, 1, 2, ..., N-1] for the row indices.
    row_indices = np.arange(num_simulations)
    # We use advanced indexing to pick the Y value from each row.
    y_samples = x_others[row_indices, second_closest_idx]

    # Step 6: Calculate E[Y], E[Y^2], and Var(Y) from the samples of Y.
    # The Law of Large Numbers states that the sample mean converges to the expectation.
    mean_y = np.mean(y_samples)
    mean_y_squared = np.mean(y_samples**2)
    
    # The variance is E[Y^2] - (E[Y])^2
    variance_y = mean_y_squared - (mean_y ** 2)
    # Alternatively, we can use np.var, which is numerically more stable:
    # variance_y = np.var(y_samples)

    # Print the result in the requested equation format
    print("\nThe variance of Y is calculated as Var(Y) = E[Y^2] - (E[Y])^2")
    print("Based on the simulation:")
    print(f"E[Y] \u2248 {mean_y}")
    print(f"E[Y^2] \u2248 {mean_y_squared}")
    print(f"Var(Y) \u2248 {mean_y_squared} - ({mean_y})^2 = {variance_y}")

# Run the calculation
find_variance_of_y()