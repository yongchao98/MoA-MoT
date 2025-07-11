import numpy as np
from fractions import Fraction

def solve_system():
    """
    Solves the given system of matrix equations and finds the sum of the first
    coordinates of the solutions.
    """

    # --- First Equation ---
    # Equation: A1*X1^2 + X1^2*B1 = C1
    # A1 is diag(5, -5), B1 is diag(6, 6), C1 is diag(-53/12, 0)
    # Let Y1 = X1^2 = [[a, b], [c, d]]
    # (5+6)a = -53/12  => 11a = -53/12
    # (5+6)b = 0       => 11b = 0
    # (-5+6)c = 0      => c = 0
    # (-5+6)d = 0      => d = 0
    # So Y1(1,1) = a = (-53/12) / 11

    # Calculate the (1,1) element of Y1 = X1^2
    y1_11_frac = Fraction(-53, 12) / 11
    y1_11 = float(y1_11_frac)

    # From X1^2 = diag(y1_11, 0), the solutions for X1 are diag(s1, 0)
    # where s1^2 = y1_11.
    # The first coordinates of solutions for X1 are s1_a and s1_b.
    s1_a = np.sqrt(complex(y1_11))
    s1_b = -s1_a

    # --- Second Equation ---
    # Equation: A2*X2^2 + X2^2*B2 = C2
    # A2 is diag(4, -5), B2 is diag(6, 6), C2 is diag(-3/11, 0)
    # Let Y2 = X2^2 = [[e, f], [g, h]]
    # (4+6)e = -3/11   => 10e = -3/11
    # (4+6)f = 0       => 10f = 0
    # (-5+6)g = 0      => g = 0
    # (-5+6)h = 0      => h = 0
    # So Y2(1,1) = e = (-3/11) / 10

    # Calculate the (1,1) element of Y2 = X2^2
    y2_11_frac = Fraction(-3, 11) / 10
    y2_11 = float(y2_11_frac)
    
    # From X2^2 = diag(y2_11, 0), the solutions for X2 are diag(s2, 0)
    # where s2^2 = y2_11.
    # The first coordinates of solutions for X2 are s2_a and s2_b.
    s2_a = np.sqrt(complex(y2_11))
    s2_b = -s2_a
    
    # The set of solutions for the first coordinate is {s1_a, s1_b, s2_a, s2_b}
    # Calculate the sum of all these first coordinates
    total_sum = s1_a + s1_b + s2_a + s2_b
    
    # Print the equation representing the sum of all first coordinates
    # The requirement is to output each number in the final equation.
    print("The individual first coordinates of the solutions are:")
    print(f"For X1: {s1_a}, {s1_b}")
    print(f"For X2: {s2_a}, {s2_b}")
    print("\nThe sum is calculated as follows:")
    print(f"({s1_a}) + ({s1_b}) + ({s2_a}) + ({s2_b}) = {total_sum}")

solve_system()
<<<0>>>