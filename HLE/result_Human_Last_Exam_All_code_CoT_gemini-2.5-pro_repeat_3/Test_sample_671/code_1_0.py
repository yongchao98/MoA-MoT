import math

def calculate_fn(n, p):
    """
    Calculates the value of the polynomial f_n(p) based on the derived formula.
    
    The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1).
    It handles the special case where p = 0.5.
    
    Args:
        n (int): The index of the polynomial in the sequence.
        p (float): The value at which to evaluate the polynomial.
    """
    # Use a small tolerance for floating point comparison
    if abs(2 * p - 1) < 1e-9:
        # Special case for p = 0.5, derived using L'Hopital's rule
        result = n * (0.5)**(n - 1)
        print(f"f_{n}({p}) is a special case (denominator is zero).")
        print(f"The value is calculated using the limit: n * (0.5)^(n-1) = {n} * {0.5**(n-1)} = {result}")
    else:
        # General formula
        p_n = p**n
        one_minus_p = 1 - p
        one_minus_p_n = one_minus_p**n
        denominator = 2 * p - 1
        numerator = p_n - one_minus_p_n
        result = numerator / denominator
        
        # Print the equation with all the numbers
        print(f"f_{n}({p}) = ({p}^{n} - ({one_minus_p})^{n}) / (2*{p} - 1)")
        print(f"        = ({p_n} - {one_minus_p_n}) / {denominator}")
        print(f"        = {numerator} / {denominator}")
        print(f"        = {result}")

if __name__ == '__main__':
    # Demonstrate the function with a few examples
    print("--- Example 1: Calculating f_5(3) ---")
    calculate_fn(5, 3)
    print("\n--- Example 2: Calculating f_15(2) ---")
    calculate_fn(15, 2)
    print("\n--- Example 3: Calculating f_4(0.5) ---")
    calculate_fn(4, 0.5)