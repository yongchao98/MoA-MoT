import math

def calculate_p_n():
    """
    Calculates P(n) using the derived closed-form formula.
    It prompts the user for n and prints the components of the formula as well as the final result.
    """
    try:
        n_str = input("Please enter a positive integer value for n: ")
        n = int(n_str)
        if n <= 0:
            print("Error: Please enter a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The closed-form formula is P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)

    # Calculate components of the formula
    two_pi = 2 * math.pi
    exponent = n * (n - 1) / 4.0
    n_factorial = float(math.factorial(n))
    sqrt_n_factorial = math.sqrt(n_factorial)

    # Calculate the final result
    try:
        numerator = math.pow(two_pi, exponent)
    except ValueError as e:
        print(f"Error calculating numerator: {e}")
        return

    if sqrt_n_factorial == 0:
        print("Error: Cannot divide by zero (sqrt(n!) is zero).")
        return

    result = numerator / sqrt_n_factorial

    # Print the explanation and the numbers in the final equation
    print(f"\nFor n = {n}, the closed-form formula for P(n) is:")
    print("P(n) = (2 * pi) ^ (n * (n - 1) / 4) / sqrt(n!)\n")
    print(f"Let's calculate the value by substituting n = {n}:")

    print(f"1. Base of the power (2 * pi): {two_pi}")
    print(f"2. Exponent (n * (n - 1) / 4): {n} * ({n} - 1) / 4 = {exponent}")
    print(f"3. Numerator ((2 * pi) ^ exponent): {numerator}")
    print(f"4. n factorial (n!): {n_factorial}")
    print(f"5. Denominator (sqrt(n!)): {sqrt_n_factorial}")

    print(f"\nFinal Result P({n}) = {numerator} / {sqrt_n_factorial} = {result}")

# Execute the function
calculate_p_n()