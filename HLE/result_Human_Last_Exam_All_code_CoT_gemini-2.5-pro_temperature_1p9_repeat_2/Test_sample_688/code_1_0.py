import math

def calculate_cn_prefactor(n):
    """
    This function calculates the system-independent prefactor c_n for the
    fully f-connected Ree-Hoover diagram contribution to the n-th virial
    coefficient, B_n.

    Args:
        n (int): The number of particles in the diagram (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Invalid input: n must be an integer greater than or equal to 2.")
        return

    numerator = n - 1
    denominator = math.factorial(n)
    value = -numerator / denominator

    # Outputting each number in the final equation as requested.
    print(f"The prefactor c_n is given by the formula: c_n = -(n-1) / n!")
    print(f"For n = {n}:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"     = -{numerator} / {denominator}")
    print(f"     = {value}")


# You can change this value to calculate the prefactor for a different n.
n_to_calculate = 4
calculate_cn_prefactor(n_to_calculate)
