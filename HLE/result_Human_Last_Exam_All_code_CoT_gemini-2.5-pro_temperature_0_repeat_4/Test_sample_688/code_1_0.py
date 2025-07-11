import math

def calculate_cn(n):
    """
    Calculates the prefactor c_n for the fully f-connected Ree-Hoover diagram
    in the virial expansion.

    The formula for c_n is: c_n = -(n-1) / n!

    Args:
        n (int): The number of nodes in the diagram (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # Calculate the numerator (n-1)
    numerator = n - 1

    # Calculate the denominator (n!)
    denominator = math.factorial(n)

    # Calculate c_n
    c_n = -numerator / denominator

    # Output the equation with all the numbers
    print(f"For n = {n}:")
    print(f"c_n = - (n - 1) / n!")
    print(f"c_{n} = - ({n} - 1) / {n}!")
    print(f"c_{n} = - {numerator} / {denominator}")
    print(f"c_{n} = {c_n}")

if __name__ == '__main__':
    # You can change this value to calculate c_n for a different n
    n_value = 4
    calculate_cn(n_value)