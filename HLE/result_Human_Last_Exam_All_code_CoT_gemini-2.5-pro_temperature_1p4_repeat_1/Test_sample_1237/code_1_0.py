import decimal

def solve_hopfield_problem():
    """
    Solves for the expression 1000*(tanh(c1)/tanh(c2) - 1)^2
    based on the given Hopfield model parameters.
    """
    # Set the precision for decimal calculations. The numbers involved are
    # extremely large (w13 ~ 10^400), so we need very high precision.
    # The result will be of the order of 10^793, so a precision of
    # 800 digits is appropriate.
    decimal.getcontext().prec = 800

    # Given model parameters as Decimal objects
    w12 = decimal.Decimal('1E5')
    w13 = decimal.Decimal('1E400')
    factor = decimal.Decimal('1000')
    one = decimal.Decimal(1)
    two = decimal.Decimal(2)

    # The decimal library does not have a built-in tanh function.
    # We calculate it using its definition: tanh(x) = (e^x - e^-x) / (e^x + e^-x)
    e_one = one.exp()
    e_neg_one = (-one).exp()
    tanh_1 = (e_one - e_neg_one) / (e_one + e_neg_one)

    # Calculate the ratio R = tanh(c1) / tanh(c2) using the derived formula:
    # R = (2 + w13 * tanh(1)) / (w12 * tanh(1))
    numerator = two + w13 * tanh_1
    denominator = w12 * tanh_1
    R = numerator / denominator

    # Calculate the final expression E = 1000 * (R - 1)^2
    term = R - one
    result = factor * (term ** 2)

    # Print the final equation with the computed values in scientific notation
    # to show each number in the final calculation.
    print(f"Based on the solvability conditions, the ratio is calculated:")
    print(f"R = tanh(c1)/tanh(c2) = (2 + w13*tanh(1)) / (w12*tanh(1)) = {R:.5e}")
    print("\nFinal equation with values:")
    print(f"{int(factor)} * ({R:.5e} - 1)^2 = {result:.5e}")

solve_hopfield_problem()