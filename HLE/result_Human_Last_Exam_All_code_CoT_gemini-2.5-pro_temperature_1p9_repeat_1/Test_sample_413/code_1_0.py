import fractions

def solve_and_print():
    """
    This script solves for the coefficients of the polynomial f(x) based on the problem's conditions
    and then computes the exact value of f(3).
    """

    # From the step-by-step derivation, we found the coefficients:
    a = fractions.Fraction(-3, 8)
    b = fractions.Fraction(-5, 8)
    c = fractions.Fraction(0)

    # The value to substitute into the function
    x_val = 3
    
    # Calculate f(3)
    val_1 = x_val**3
    val_2 = a * (x_val**2)
    val_3 = b * x_val
    
    result = val_1 + val_2 + val_3 + c

    # Print the calculation steps for clarity
    print(f"Based on the analysis, the polynomial is f(x) = x^3 + ({a})x^2 + ({b})x.")
    print("Now we compute f(3):")
    print(f"f(3) = ({x_val})^3 + ({a}) * ({x_val})^2 + ({b}) * ({x_val})")
    print(f"     = {val_1} + ({a}) * 9 + ({b}) * 3")
    print(f"     = {val_1} + ({val_2}) + ({val_3})")
    print(f"     = 27 - 27/8 - 15/8")
    print(f"     = 27 - 42/8")
    print(f"     = 27 - 21/4")
    print(f"     = 108/4 - 21/4")
    print(f"     = {result}")

solve_and_print()
<<<87/4>>>