import math

def solve_for_c(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    Args:
        d (int): The dimension d, must be >= 1.
        p (float): The parameter p, must be > 0.

    Returns:
        float: The smallest constant c.
    """
    if d < 1 or p <= 0:
        raise ValueError("d must be >= 1 and p must be > 0")

    # The exponent is max(0, 1/2 - 1/p)
    exponent = max(0, 0.5 - 1/p)
    
    # The constant c is d raised to the power of the exponent
    c = d ** exponent
    
    return c

if __name__ == '__main__':
    # Example values for d and p
    d_val = 10
    p_val = 4.0

    # Calculate c
    c_val = solve_for_c(d_val, p_val)

    # Output the equation with the calculated numbers
    print(f"For d = {d_val} and p = {p_val}:")
    exponent_val = max(0, 0.5 - 1/p_val)
    print(f"c = d^max(0, 1/2 - 1/p) = {d_val}^max(0, 0.5 - {1/p_val:.3f}) = {d_val}^{exponent_val:.3f} = {c_val:.4f}")
    
    # Another example for p < 2
    d_val_2 = 10
    p_val_2 = 1.5

    c_val_2 = solve_for_c(d_val_2, p_val_2)
    print("\n" + "="*20 + "\n")
    print(f"For d = {d_val_2} and p = {p_val_2}:")
    exponent_val_2 = max(0, 0.5 - 1/p_val_2)
    print(f"c = d^max(0, 1/2 - 1/p) = {d_val_2}^max(0, 0.5 - {1/p_val_2:.3f}) = {d_val_2}^{exponent_val_2:.3f} = {c_val_2:.4f}")
