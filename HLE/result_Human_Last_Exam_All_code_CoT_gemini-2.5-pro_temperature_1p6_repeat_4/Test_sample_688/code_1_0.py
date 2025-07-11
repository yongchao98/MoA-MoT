import math

def calculate_cn_prefactor(n):
    """
    Calculates and prints the system-independent prefactor c_n for the fully
    f-connected Ree-Hoover diagram with n nodes.

    The derived formula for the prefactor is c_n = -(n-1) / n!

    Args:
        n (int): The number of nodes (particles) in the diagram.
                 Must be an integer greater than or equal to 2.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # Calculate the numerator of the fraction
    numerator = n - 1

    # Calculate the denominator of the fraction (n!)
    denominator = math.factorial(n)

    # Calculate the final result
    result = -numerator / denominator

    # Print the equation and the result, showing each number involved.
    print(f"For n = {n}:")
    print(f"The formula is c_n = - (n - 1) / n!")
    print(f"c_{n} = - ({n} - 1) / {n}!")
    print(f"c_{n} = -{numerator}/{denominator}")
    print(f"c_{n} = {result}")
    print("-" * 20)

if __name__ == '__main__':
    # Demonstrate the calculation for a few values of n
    print("Calculating the prefactor c_n for several values of n.")
    print("-" * 20)
    for i in range(2, 6):
        calculate_cn_prefactor(i)
