import numpy as np

def calculate_variance_of_y():
    """
    This function calculates the variance of Y using a Monte Carlo simulation.

    Let X_1, X_2, X_3, and X_4 be i.i.d. random variables uniformly distributed
    on the interval [0, 1]. Y is the second closest value to X_1 among
    the sample points X_2, X_3, and X_4.

    The variance is determined by generating a large number of samples and
    computing their sample variance.
    """
    # Set a seed for reproducibility of the random numbers
    np.random.seed(42)

    # Number of simulations to run. A larger number yields a more accurate estimate.
    num_simulations = 5_000_000

    # List to store the simulated values of Y
    y_samples = []

    # Main simulation loop
    for _ in range(num_simulations):
        # Step a: Generate four random values from U[0, 1]
        x = np.random.rand(4)
        x1 = x[0]
        other_xs = x[1:]  # This array holds X_2, X_3, X_4

        # Step b: Calculate the absolute distances from x1
        distances = np.abs(other_xs - x1)

        # Step c: Find the value of X_i corresponding to the second smallest distance.
        # We can use argsort to find the indices that would sort the distances array.
        # The index of the second smallest distance is at position 1 of the sorted indices.
        second_closest_index = np.argsort(distances)[1]
        y_sample = other_xs[second_closest_index]

        # Step d: Store the sample of Y
        y_samples.append(y_sample)

    # Step e: Calculate the variance of the collected Y samples.
    # We also calculate the mean E[Y] and the mean of squares E[Y^2].
    mean_y = np.mean(y_samples)
    mean_y_squared = np.mean(np.square(y_samples))
    variance_y = np.var(y_samples) # Uses N in denominator, consistent with E[Y^2] - E[Y]^2

    print(f"Simulation based on {num_simulations:,} samples.")
    print("-" * 40)
    print("The variance of Y is given by the formula: Var(Y) = E[Y^2] - (E[Y])^2")
    print("\nEstimated values from the simulation:")
    print(f"E[Y] ≈ {mean_y:.6f}")
    print(f"E[Y^2] ≈ {mean_y_squared:.6f}")
    
    print("\nPlugging these values into the formula:")
    print(f"Var(Y) ≈ {mean_y_squared:.6f} - ({mean_y:.6f})^2")
    print(f"Var(Y) ≈ {mean_y_squared:.6f} - {mean_y**2:.6f}")
    print(f"Var(Y) ≈ {variance_y:.6f}")
    print("-" * 40)
    print(f"The final estimated variance of Y is: {variance_y:.6f}")

# Execute the calculation
calculate_variance_of_y()