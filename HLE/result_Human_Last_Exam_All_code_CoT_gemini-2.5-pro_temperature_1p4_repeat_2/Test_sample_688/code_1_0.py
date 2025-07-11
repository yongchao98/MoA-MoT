import math

def calculate_cn(n):
    """
    Calculates the system-independent prefactor c_n for a given integer n >= 2.

    The formula for c_n is derived from the Mayer cluster expansion for the virial
    coefficients and is given by c_n = -1 / (n * (n-2)!).

    Args:
        n (int): The number of nodes in the diagram, must be >= 2.
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer greater than or equal to 2. Received n={n}.")
        return

    try:
        # Calculate the components of the formula
        n_minus_2_fact = math.factorial(n - 2)
        denominator = n * n_minus_2_fact
        cn_value = -1 / denominator

        # Print the detailed calculation as requested
        print(f"Calculation for n = {n}:")
        print(f"c_n = -1 / (n * (n-2)!)")
        print(f"c_{n} = -1 / ({n} * ({n}-2)!)")
        print(f"c_{n} = -1 / ({n} * {n-2}!)")
        print(f"c_{n} = -1 / ({n} * {n_minus_2_fact})")
        print(f"c_{n} = -1 / {denominator}")
        print(f"c_{n} = {cn_value:.10f}")
        print("-" * 30)

    except ValueError:
        print(f"Error calculating factorial for n={n}. Input must be a non-negative integer.")

if __name__ == '__main__':
    # Demonstrate the calculation for a few values of n
    calculate_cn(2)
    calculate_cn(3)
    calculate_cn(4)
    calculate_cn(5)
    calculate_cn(6)