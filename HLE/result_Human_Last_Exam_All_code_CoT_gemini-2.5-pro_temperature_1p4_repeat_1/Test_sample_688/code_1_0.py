import math

def calculate_cn_prefactor(n: int):
    """
    Calculates the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram contribution to the n-th virial coefficient.

    The formula for the prefactor is c_n = -(n-1) / n!

    Args:
        n: The order of the virial coefficient (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # Calculate the components of the formula
    numerator = n - 1
    denominator = math.factorial(n)
    cn_value = -numerator / denominator

    # Print the result, showing the numbers in the equation
    print(f"For n = {n}:")
    print(f"The formula is: c_n = - (n - 1) / n!")
    # Print the equation with numbers substituted
    print(f"c_{n} = - ({n} - 1) / {n}!")
    print(f"c_{n} = - {numerator} / {denominator}")
    # Print the final numerical value
    print(f"c_{n} = {cn_value}\n")

if __name__ == '__main__':
    # Demonstrate the calculation for a few values of n
    # The prefactor is defined for n >= 2
    calculate_cn_prefactor(2)
    calculate_cn_prefactor(3)
    calculate_cn_prefactor(4)
    calculate_cn_prefactor(5)
