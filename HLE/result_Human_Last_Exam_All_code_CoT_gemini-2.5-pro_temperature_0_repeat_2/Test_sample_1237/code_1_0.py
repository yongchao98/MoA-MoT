import decimal

def solve_hopfield_problem():
    """
    Calculates the value of the expression 1000 * (tanh(c1)/tanh(c2) - 1)^2
    based on the provided Hopfield model parameters.
    """
    # Set precision for decimal calculations. The numbers involved are very large,
    # so we need a high precision to handle them. 800 digits should be sufficient.
    decimal.getcontext().prec = 800

    # Define the given constants as Decimal objects for high-precision arithmetic
    w12 = decimal.Decimal('1e5')
    w13 = decimal.Decimal('1e400')
    val_2 = decimal.Decimal('2')
    val_1000 = decimal.Decimal('1000')

    # The decimal module does not have a built-in tanh function.
    # We can define it using tanh(x) = (exp(x) - exp(-x)) / (exp(x) + exp(-x)).
    one = decimal.Decimal(1)
    exp1 = one.exp()
    exp_minus_1 = (-one).exp()
    tanh1 = (exp1 - exp_minus_1) / (exp1 + exp_minus_1)

    # Calculate the ratio R = tanh(c1)/tanh(c2)
    # R = (2 + w13 * tanh1) / (w12 * tanh1)
    ratio_numerator = val_2 + w13 * tanh1
    ratio_denominator = w12 * tanh1
    ratio = ratio_numerator / ratio_denominator

    # Calculate the final expression: 1000 * (R - 1)^2
    final_result = val_1000 * (ratio - 1)**2

    # Print the final equation with all the numbers substituted, as requested.
    # We use scientific notation for clarity due to the large numbers.
    print("The expression to calculate is: 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2")
    print(f"Derived ratio: tanh(c1)/tanh(c2) = (2 + w13*tanh(1)) / (w12*tanh(1))")
    print("Substituting the given values:")
    print(f"1000 * ( (2 + {w13.to_eng_string()} * {tanh1}) / ({w12.to_eng_string()} * {tanh1}) - 1 )^2 = {final_result.to_eng_string()}")

solve_hopfield_problem()