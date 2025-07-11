import math

def calculate_lewis_weight_constant(d, p):
    """
    Calculates the smallest constant c such that for any matrix A,
    ||W^(1/2-1/p)Ax||_2 <= c ||Ax||_p for every x.

    Args:
        d (int): The dimension of the column space of matrix A.
        p (float): The L_p norm parameter.

    Returns:
        float: The value of the constant c.
    """
    if d <= 0 or not isinstance(d, int):
        raise ValueError("Dimension 'd' must be a positive integer.")
    if p <= 0:
        raise ValueError("Norm parameter 'p' must be a positive number.")

    # The formula for the constant c is d**|1/2 - 1/p|
    exponent_val = 0.5 - (1 / p)
    abs_exponent_val = abs(exponent_val)
    c = d ** abs_exponent_val

    # Outputting the numbers in the final equation as requested
    print(f"For d = {d} and p = {p}:")
    print("The formula for the constant c is: c = d**|1/2 - 1/p|")
    print(f"c = {d}**|1/2 - 1/{p}|")
    print(f"c = {d}**|{0.5 - 1/p}|")
    print(f"c = {d}**{abs_exponent_val}")
    print(f"The final value is: c = {c}")
    return c

if __name__ == '__main__':
    # Example 1: d=4, p=4
    print("--- Example 1 ---")
    d1 = 4
    p1 = 4.0
    calculate_lewis_weight_constant(d1, p1)

    print("\n" + "--- Example 2 ---")
    # Example 2: d=9, p=1
    d2 = 9
    p2 = 1.0
    calculate_lewis_weight_constant(d2, p2)

    print("\n" + "--- Example 3 ---")
    # Example 3: d=1, p=2 (should be 1)
    d3 = 1
    p3 = 2.0
    calculate_lewis_weight_constant(d3, p3)
