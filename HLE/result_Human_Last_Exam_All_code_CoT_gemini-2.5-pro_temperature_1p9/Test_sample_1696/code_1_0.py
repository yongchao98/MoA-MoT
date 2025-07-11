import math

def calculate_P_n_closed_form(n):
    """
    Calculates the closed-form formula for P(n) and prints the result.
    The closed-form for P(n) is 1/sqrt(n!).
    
    Args:
    n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    try:
        n_factorial = math.factorial(n)
        result = 1 / math.sqrt(n_factorial)
        
        print(f"For n = {n}:")
        print(f"P({n}) = 1 / sqrt({n}!)")
        print(f"P({n}) = 1 / sqrt({n_factorial})")
        print(f"P({n}) = {result}")

    except ValueError:
        print(f"Factorial is not defined for negative numbers.")
    except OverflowError:
        print(f"Cannot calculate factorial for n={n}. The number is too large.")

# Example usage of the function for n=5
n_value = 5
calculate_P_n_closed_form(n_value)

# Example usage of the function for n=10
n_value_2 = 10
calculate_P_n_closed_form(n_value_2)
