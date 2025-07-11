import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
    """
    # The problem reduces to finding the limit L = 2 / ln(alpha),
    # where alpha is the larger root of the characteristic equation r^2 - 13r + 1 = 0.
    # The equation for alpha is derived from the recurrence relation that generates the solutions.

    # Step 1: Define the equation to be solved for the final value.
    # The value is floor(10000 * (2 / ln( (13 + sqrt(165)) / 2 )))
    print("The expression to calculate is: floor(10000 * (2 / ln((13 + sqrt(165)) / 2)))")
    print("-" * 20)

    # Step 2: Calculate the components of the equation.
    # Calculate sqrt(165)
    val_sqrt_165 = math.sqrt(165)
    print(f"First, we calculate the value of sqrt(165):")
    print(f"sqrt(165) = {val_sqrt_165}")
    print("-" * 20)

    # Calculate alpha = (13 + sqrt(165)) / 2
    val_alpha = (13 + val_sqrt_165) / 2
    print(f"Next, we calculate the constant alpha = (13 + sqrt(165)) / 2:")
    print(f"alpha = (13 + {val_sqrt_165}) / 2 = {val_alpha}")
    print("-" * 20)

    # Calculate ln(alpha)
    val_ln_alpha = math.log(val_alpha)
    print(f"Then, we take the natural logarithm of alpha:")
    print(f"ln(alpha) = ln({val_alpha}) = {val_ln_alpha}")
    print("-" * 20)

    # Calculate the limit L = 2 / ln(alpha)
    limit_L = 2 / val_ln_alpha
    print(f"The limit value is L = 2 / ln(alpha):")
    print(f"L = 2 / {val_ln_alpha} = {limit_L}")
    print("-" * 20)

    # Step 3: Calculate the final result.
    final_value = 10000 * limit_L
    print(f"Finally, we multiply the limit by 10^4:")
    print(f"10000 * L = 10000 * {limit_L} = {final_value}")
    print("-" * 20)

    # Step 4: Find the integer part.
    integer_part = math.floor(final_value)
    print(f"The integer part of the result is:")
    print(f"floor({final_value}) = {integer_part}")

solve()