import math

def calculate_cn_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    Ree-Hoover diagram in the virial expansion.

    The formula for the prefactor is c_n = -(n-1) / n!

    Args:
        n (int): The number of nodes in the diagram (must be >= 2).

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return None

    # Calculate the numerator of the fraction
    numerator = n - 1

    # Calculate the denominator of the fraction (n!)
    denominator = math.factorial(n)

    # Calculate the prefactor c_n
    c_n = -numerator / denominator

    # Print the result, showing the components of the calculation
    print(f"The system-independent prefactor c_n is given by the formula: c_n = -(n-1)/n!")
    print(f"For n = {n}:")
    print(f"c_{n} = -({n} - 1) / {n}!")
    print(f"c_{n} = -{numerator} / {denominator}")
    print(f"c_{n} = {c_n}")
    
    return c_n

if __name__ == '__main__':
    # Example: Calculate the prefactor for n=4
    n_value = 4
    calculate_cn_prefactor(n_value)
    
    print("\n" + "="*20 + "\n")

    # Example: Calculate the prefactor for n=5
    n_value = 5
    calculate_cn_prefactor(n_value)