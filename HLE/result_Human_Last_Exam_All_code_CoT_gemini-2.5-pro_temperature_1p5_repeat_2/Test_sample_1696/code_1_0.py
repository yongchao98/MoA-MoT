import math

def calculate_P(n):
    """
    Calculates the closed-form formula for P(n) and prints the steps.
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    print(f"For n = {n}:")

    # Calculate the exponent part
    exponent_numerator = n * (n - 1)
    exponent_denominator = 4
    exponent_val = exponent_numerator / exponent_denominator
    print(f"The exponent is n*(n-1)/4 = {n}*({n}-1)/{exponent_denominator} = {exponent_numerator}/{exponent_denominator} = {exponent_val}")
    
    # Calculate n!
    try:
        n_factorial = math.factorial(n)
        print(f"The factorial is n! = {n}! = {n_factorial}")
    except ValueError:
         print(f"Error: Factorial not defined for n={n}")
         return


    # Assemble the formula string
    formula_str = f"P({n}) = (2*pi)^({exponent_val}) / sqrt({n_factorial})"
    print("The formula is: " + formula_str)

    # Calculate the final numerical value
    try:
        numerator = (2 * math.pi)**(exponent_val)
        denominator = math.sqrt(n_factorial)
        result = numerator / denominator
        print(f"The numerical value is P({n}) = {result}")
    except OverflowError:
        print(f"Numerical result for P({n}) is too large to represent as a float.")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")

# Example usage with n=10, as the user did not specify a value for n.
# You can change the value of n to any positive integer.
n = 10
calculate_P(n)