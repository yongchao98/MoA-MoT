import numpy as np

def solve_model():
    """
    This function devises a maximally parsimonious model for the given data,
    estimates its parameters, and prints the final equation.
    """
    # Step 1: Define the dataset
    x = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Step 2: Propose and fit the model.
    # A scatter plot of the data suggests a relationship of the form y = sqrt(f(x)).
    # The simplest parsimonious model of this type is y = sqrt(a*x + b).
    # We can linearize this to y^2 = a*x + b and find the parameters 'a' and 'b'
    # using a standard linear least-squares fit.
    y_squared = y**2

    # Step 3: Use numpy.polyfit to perform a linear regression on y^2 vs x.
    # This finds the coefficients [a, b] that minimize the sum of squared residuals.
    params = np.polyfit(x, y_squared, 1)
    a, b = params

    # Step 4: Report the results.
    # Format the parameters to 3 significant digits.
    a_str = f"{a:.3g}"
    b_val = b
    
    # Format the final equation string with signs handled correctly.
    if b_val < 0:
        sign = "-"
        b_abs_str = f"{abs(b_val):.3g}"
        equation = f"y = sqrt({a_str}*x {sign} {b_abs_str})"
    else:
        sign = "+"
        b_abs_str = f"{b_val:.3g}"
        equation = f"y = sqrt({a_str}*x {sign} {b_abs_str})"

    print("The maximally parsimonious model that optimally predicts y is:")
    print(equation)
    print("\nThe estimated parameter values are:")
    print(f"a = {a_str}")
    print(f"b = {b_val:.3g}")

solve_model()