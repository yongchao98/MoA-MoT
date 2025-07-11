import cmath
from fractions import Fraction

def solve_and_sum():
    """
    Solves the two matrix equations to find the first coordinates of the solutions
    and calculates their sum.
    """
    
    # First equation
    # (A1 + 6I) * Y1 = C1
    # Y1 = X1^2
    c1_11 = Fraction(-53, 12)
    a1_11 = 5
    b_11 = 6
    # (a1_11 + b_11) * y1_11 = c1_11
    y1_11 = c1_11 / (a1_11 + b_11)
    
    # The first coordinate of X1 is the square root of y1_11
    # x1_11^2 = y1_11
    x1_sol1 = cmath.sqrt(y1_11)
    x1_sol2 = -x1_sol1
    
    print("For the first equation:")
    print(f"The (1,1) element of X1^2 is {y1_11.numerator}/{y1_11.denominator}")
    print(f"The first coordinates of the solutions for X1 are: {x1_sol1} and {x1_sol2}")
    print("-" * 20)

    # Second equation
    # (A2 + 6I) * Y2 = C2
    # Y2 = X2^2
    c2_11 = Fraction(-3, 11)
    a2_11 = 4
    # b_11 is the same (6)
    # (a2_11 + b_11) * y2_11 = c2_11
    y2_11 = c2_11 / (a2_11 + b_11)
    
    # The first coordinate of X2 is the square root of y2_11
    # x2_11^2 = y2_11
    x2_sol1 = cmath.sqrt(y2_11)
    x2_sol2 = -x2_sol1

    print("For the second equation:")
    print(f"The (1,1) element of X2^2 is {y2_11.numerator}/{y2_11.denominator}")
    print(f"The first coordinates of the solutions for X2 are: {x2_sol1} and {x2_sol2}")
    print("-" * 20)
    
    # Sum of all first coordinates
    total_sum = x1_sol1 + x1_sol2 + x2_sol1 + x2_sol2
    
    # Output the final equation with each number
    print("The sum of the first coordinates of all solutions is:")
    print(f"({x1_sol1}) + ({x1_sol2}) + ({x2_sol1}) + ({x2_sol2}) = {total_sum}")

solve_and_sum()