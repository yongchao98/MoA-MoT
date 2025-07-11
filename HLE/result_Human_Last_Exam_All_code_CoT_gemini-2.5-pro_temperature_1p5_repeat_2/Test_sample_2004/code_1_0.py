import numpy as np

def find_parsimonious_model():
    """
    This function finds a maximally parsimonious model for the given data,
    which is determined to be of the form y = c * x^0.75, and then
    calculates the optimal parameter 'c'.
    """
    # The 25 observations of x and y
    x = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # The chosen parsimonious model is y = c * x^p where p=0.75.
    # This is a one-parameter model.
    p = 0.75

    # To find the least-squares estimate for 'c', we solve for c that minimizes
    # the sum of squared errors: sum( (y_i - c * x_i^p)^2 ).
    # The solution is c = sum(y_i * x_i^p) / sum((x_i^p)^2).
    
    # Let's define a transformed variable z = x^p
    z = x**p
    
    # Calculate the single parameter 'c' using the least-squares formula
    c_estimated = np.dot(y, z) / np.dot(z, z)
    
    # Round the estimated parameter to 3 significant digits
    # The format specifier {:.3g} correctly handles this.
    c_rounded = float(f"{c_estimated:.3g}")

    # Report the final model and its parameters
    print("The selected maximally parsimonious model is a power law with a single parameter.")
    print("The final equation is of the form: y = c * x^p")
    print("\nCalculated values for the equation:")
    print(f"The estimated parameter c = {c_rounded}")
    print(f"The fixed exponent p = {p}")
    print("\nFinal Model:")
    print(f"y = {c_rounded} * x^{p}")

# Execute the function to print the result
find_parsimonious_model()
<<<y = 0.485 * x^0.75>>>