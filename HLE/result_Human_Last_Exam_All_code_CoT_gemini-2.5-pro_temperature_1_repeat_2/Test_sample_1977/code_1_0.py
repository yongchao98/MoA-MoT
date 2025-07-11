import math

def calculate_and_print_norm(n):
    """
    Calculates and prints the 1-norm of the correlation matrix T for a given even integer n.
    The function first checks if n is even. It then computes the norm using the derived formula:
    ||T||_1 = (2**(n+1) * 3**n + 3**(n+1) - 2**(n+1) - 1) / (1 + 3**n)
    It prints the calculation steps, showing how the numbers are plugged into the formula.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print(f"Error: The problem is defined for even non-negative integers. Your input was n = {n}.")
        return

    # To handle potentially large numbers, all calculations are done with Python's arbitrary-precision integers.
    
    # Numerator of the formula: 2**(n+1) * 3**n + 3**(n+1) - 2**(n+1) - 1
    term1_num = pow(2, n + 1) * pow(3, n)
    term2_num = pow(3, n + 1)
    term3_num = pow(2, n + 1)
    term4_num = 1
    numerator = term1_num + term2_num - term3_num - term4_num

    # Denominator of the formula: 1 + 3**n
    denominator = 1 + pow(3, n)

    # The final result. It might be a fraction, so we perform floating-point division at the end.
    result = numerator / denominator

    print(f"For n = {n}:")
    print("The formula for the 1-norm of the correlation matrix T is:")
    print(f"||T||_1 = (2^({n}+1) * 3^{n} + 3^({n}+1) - 2^({n}+1) - 1) / (1 + 3^{n})")
    print()
    print("Plugging in the value of n:")
    print(f"        = ({pow(2, n+1)} * {pow(3, n)} + {pow(3, n+1)} - {pow(2, n+1)} - 1) / (1 + {pow(3, n)})")
    print()
    print("Evaluating the terms:")
    print(f"        = ({term1_num} + {term2_num} - {term3_num} - {term4_num}) / {denominator}")
    print()
    print("Calculating the final fraction:")
    print(f"        = {numerator} / {denominator}")
    print()
    print("The final result is:")
    print(f"        = {result}")
    print("-" * 30)

if __name__ == "__main__":
    # You can test the function for any even integer n.
    # Here are some examples:
    calculate_and_print_norm(2)
    calculate_and_print_norm(4)
    calculate_and_print_norm(10)