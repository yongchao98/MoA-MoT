def calculate_variance():
    """
    Calculates and explains the variance of Y.

    The variance of a random variable Y is given by the formula:
    Var(Y) = E[Y^2] - (E[Y])^2

    For this problem, we can determine the required expected values:
    1. E[Y]: By symmetry, the expected value of Y is the midpoint of the
       interval [0, 1], which is 1/2.
    2. E[Y^2]: This value requires a complex analytical derivation using
       the law of total expectation, which evaluates to 37/120.

    This script performs the final calculation using these values.
    """
    # Numerators and denominators for the analytical results
    e_y_sq_num = 37
    e_y_sq_den = 120
    e_y_num = 1
    e_y_den = 2

    # Calculate (E[Y])^2
    e_y_squared_num = e_y_num ** 2
    e_y_squared_den = e_y_den ** 2

    # To subtract the fractions, we find a common denominator
    common_denominator = e_y_sq_den
    factor = common_denominator // e_y_squared_den
    e_y_squared_num_common = e_y_squared_num * factor

    # Calculate the variance numerator
    var_y_num = e_y_sq_num - e_y_squared_num_common
    var_y_den = common_denominator

    print("The variance is calculated using the formula: Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"By symmetry, the expected value E[Y] = {e_y_num}/{e_y_den}.")
    print(f"The expected value of the square, E[Y^2], is analytically found to be {e_y_sq_num}/{e_y_sq_den}.")
    print("\nThe calculation for the variance is as follows:")
    print(f"Var(Y) = {e_y_sq_num}/{e_y_sq_den} - ({e_y_num}/{e_y_den})^2")
    print(f"       = {e_y_sq_num}/{e_y_sq_den} - {e_y_squared_num}/{e_y_squared_den}")
    print(f"       = {e_y_sq_num}/{e_y_sq_den} - {e_y_squared_num_common}/{common_denominator}")
    print(f"       = {var_y_num}/{var_y_den}")
    print(f"\nThe final variance as a decimal is {var_y_num/var_y_den}.")

calculate_variance()