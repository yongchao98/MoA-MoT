import math

def solve_for_c(d, p):
    """
    Calculates the smallest constant c for the given inequality.

    Args:
        d (int): The dimension d, must be a positive integer.
        p (float): The p-norm parameter, must be positive.

    Returns:
        float: The value of the constant c.
    """
    if d <= 0 or not isinstance(d, int):
        raise ValueError("d must be a positive integer.")
    if p <= 0:
        raise ValueError("p must be a positive number.")

    exponent = 0.5 - 1.0 / p
    
    # The formula is c = max(1, d^exponent)
    if exponent > 0:
        c = d ** exponent
    else:
        c = 1.0
        
    return c

def main():
    """
    Main function to demonstrate the calculation for example values of d and p.
    """
    # The problem fixes d. We choose an example value.
    d = 10
    
    # We choose an example value for p > 2.
    p1 = 4.0
    
    c1 = solve_for_c(d, p1)
    
    # Output the explanation and the result for the first case
    print("The smallest constant c is given by the formula: c = max(1, d^(1/2 - 1/p))")
    print("\n--- Example 1 (p > 2) ---")
    print(f"For d = {d} and p = {p1}:")
    
    exponent1 = 0.5 - 1.0 / p1
    power_val1 = d ** exponent1
    
    print(f"The equation for c is: c = max(1, {d}^(1/2 - 1/{p1}))")
    print(f"c = max(1, {d}^({exponent1}))")
    print(f"c = max(1, {power_val1:.4f})")
    print(f"The smallest value for c is: {c1:.4f}")

    # We choose another example value for p <= 2.
    p2 = 1.5
    
    c2 = solve_for_c(d, p2)
    
    # Output the explanation and the result for the second case
    print("\n--- Example 2 (p <= 2) ---")
    print(f"For d = {d} and p = {p2}:")
    
    exponent2 = 0.5 - 1.0 / p2
    power_val2 = d ** exponent2
    
    print(f"The equation for c is: c = max(1, {d}^(1/2 - 1/{p2}))")
    print(f"c = max(1, {d}^({exponent2:.4f}))")
    print(f"c = max(1, {power_val2:.4f})")
    print(f"The smallest value for c is: {c2:.4f}")


if __name__ == "__main__":
    main()