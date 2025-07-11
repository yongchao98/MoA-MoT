import math

def calculate_lewis_weights_inequality_constant(d, p):
    """
    Calculates the smallest constant c for the Lewis weights inequality.

    The inequality is ||W^(1/2-1/p)Ax||_2 <= c ||Ax||_p, where d is the number
    of columns of A. The constant c is given by the formula d^|1/2 - 1/p|.

    Args:
        d (int): The dimension d, a positive integer.
        p (float): The parameter p, a positive real number.

    Returns:
        float: The value of the constant c.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(p, (int, float)) or p <= 0:
        raise ValueError("p must be a positive real number.")

    exponent = abs(0.5 - 1/p)
    c = d ** exponent
    
    # The prompt asks to output each number in the final equation.
    print(f"For d = {d} and p = {p}:")
    print(f"The smallest constant c is given by the formula: c = d^|1/2 - 1/p|")
    print(f"c = {d}^|0.5 - 1/{p}|")
    print(f"c = {d}^|0.5 - {1/p}|")
    print(f"c = {d}^{exponent}")
    print(f"c = {c}")
    return c

if __name__ == '__main__':
    # You can change these values to test with different d and p
    # d must be a positive integer
    d_val = 10
    # p must be a positive float
    p_val = 4.0
    
    # Let's test with another case, p < 2
    # d_val = 4
    # p_val = 1.0

    calculate_lewis_weights_inequality_constant(d_val, p_val)
