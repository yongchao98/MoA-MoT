import math

def compute_and_print_cn(n):
    """
    Computes and prints the prefactor c_n for a given n.

    The prefactor c_n is part of the virial expansion of the equation of state.
    It is given by the formula: c_n = -(n-1)/n!

    Args:
        n (int): The number of nodes in the diagram (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2. Received: {n}")
        return

    # Calculate numerator and denominator
    numerator = n - 1
    denominator = math.factorial(n)

    # To display a simplified fraction, find the greatest common divisor
    common_divisor = math.gcd(numerator, denominator)
    simple_numerator = numerator // common_divisor
    simple_denominator = denominator // common_divisor

    # The float value of the prefactor
    value = -numerator / denominator

    # Print the results in a clear, step-by-step format
    print(f"For n = {n}:")
    print(f"The formula for the prefactor is c_n = -(n - 1) / n!")
    print(f"Substituting n = {n}, we get:")
    # Outputting each number in the final equation
    print(f"c_{n} = -({n} - 1) / {n}! = -{numerator} / {denominator}")

    # Show simplified fraction if different
    if common_divisor > 1:
        print(f"Simplified fraction: c_{n} = -{simple_numerator} / {simple_denominator}")

    # Show decimal value
    print(f"Decimal value: c_{n} = {value}")
    print("-" * 25)

if __name__ == '__main__':
    # Calculate and print the prefactor for a few example values of n
    compute_and_print_cn(2)
    compute_and_print_cn(3)
    compute_and_print_cn(4)
    compute_and_print_cn(5)
    compute_and_print_cn(6)