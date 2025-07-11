import math
from decimal import Decimal, getcontext

def solve_hopfield_problem():
    """
    Calculates the value of the expression 1000 * (tanh(c1)/tanh(c2) - 1)^2
    based on the provided Hopfield model equations and parameters.
    """
    # The magnitude of the result is approximately 10^793.
    # Set precision high enough to represent all digits of the result.
    getcontext().prec = 800

    # Given parameters as Decimal objects for high-precision arithmetic
    w13 = Decimal(10)**400
    w12 = Decimal(10)**5
    
    # Constants
    const_2 = Decimal(2)
    const_1 = Decimal(1)
    const_1000 = Decimal(1000)

    # Calculate tanh(1) using Decimal precision
    tanh1 = Decimal(math.tanh(1))

    # Calculate the ratio tanh(c1)/tanh(c2)
    # ratio = (2 + w13 * tanh(1)) / (w12 * tanh(1))
    ratio_numerator = const_2 + w13 * tanh1
    ratio_denominator = w12 * tanh1
    ratio = ratio_numerator / ratio_denominator

    # Calculate the term inside the square: (ratio - 1)
    term_in_square = ratio - const_1

    # Calculate the final result
    result = const_1000 * (term_in_square ** 2)

    # Per instructions, output each number in the final equation.
    # The equation is: 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2 = result
    # We will print the values of the components leading to the final result.
    print(f"1000 * ( ( (2 + {w13} * {tanh1}) / ({w12} * {tanh1}) ) - 1 )^2 = {result}")

solve_hopfield_problem()
