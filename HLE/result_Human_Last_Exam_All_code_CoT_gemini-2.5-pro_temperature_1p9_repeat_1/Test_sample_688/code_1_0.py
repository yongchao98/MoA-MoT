import math

def calculate_cn(n):
    """
    Calculates the system-independent prefactor c_n for a given n.

    The formula for c_n is:
    c_n = -(n - 1) / n!

    Args:
        n (int): The number of nodes (particles), must be >= 2.

    Returns:
        float: The value of the prefactor c_n.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")

    numerator = -(n - 1)
    denominator = math.factorial(n)
    cn = numerator / denominator
    
    print(f"For n = {n}:")
    print(f"The formula for the prefactor is c_n = -(n - 1) / n!")
    print(f"Numerator (-(n - 1)): -({n} - 1) = {numerator}")
    print(f"Denominator (n!): {n}! = {denominator}")
    print(f"c_{n} = {numerator} / {denominator} = {cn}")
    return cn

if __name__ == '__main__':
    # Example: Calculate c_n for n=4
    calculate_cn(4)
    print("\n----------------\n")
    # Example: Calculate c_n for n=5
    calculate_cn(5)