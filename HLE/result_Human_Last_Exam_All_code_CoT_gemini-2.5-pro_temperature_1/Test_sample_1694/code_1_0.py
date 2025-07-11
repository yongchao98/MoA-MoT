import cmath
from fractions import Fraction

def solve_and_sum():
    """
    Solves the system of matrix equations and calculates the sum of the first coordinates of the solutions.
    """
    # --- Equation 1 ---
    a1_11 = 5
    a1_22 = -5
    b_11 = 6
    b_22 = 6
    c1_11 = Fraction(-53, 12)
    c1_22 = 0

    # Solve for Y1 = X1^2
    # Since all matrices are diagonal, Y1 is diagonal.
    # (a1_11 + b_11) * y1_11 = c1_11
    # (a1_22 + b_22) * y1_22 = c1_22
    y1_11_coeff = a1_11 + b_11
    y1_11 = c1_11 / y1_11_coeff
    y1_22_coeff = a1_22 + b_22
    y1_22 = c1_22 / y1_22_coeff

    print("--- Solving for X1 ---")
    print(f"The equation for the (1,1) element of X1^2 is: ({a1_11} + {b_11}) * (X1^2)_11 = {c1_11}")
    print(f"This gives (X1^2)_11 = {c1_11} / {y1_11_coeff} = {y1_11}")
    
    # Find the first coordinates of X1
    x1_11_sol1 = cmath.sqrt(y1_11)
    x1_11_sol2 = -x1_11_sol1
    print("The two possible values for the first coordinate of X1, (X1)_11, are:")
    print(f"  1. {x1_11_sol1}")
    print(f"  2. {x1_11_sol2}\n")

    # --- Equation 2 ---
    a2_11 = 4
    a2_22 = -5
    c2_11 = Fraction(-3, 11)
    c2_22 = 0

    # Solve for Y2 = X2^2
    y2_11_coeff = a2_11 + b_11
    y2_11 = c2_11 / y2_11_coeff
    y2_22_coeff = a2_22 + b_22
    y2_22 = c2_22 / y2_22_coeff
    
    print("--- Solving for X2 ---")
    print(f"The equation for the (1,1) element of X2^2 is: ({a2_11} + {b_11}) * (X2^2)_11 = {c2_11}")
    print(f"This gives (X2^2)_11 = {c2_11} / {y2_11_coeff} = {y2_11}")

    # Find the first coordinates of X2
    x2_11_sol1 = cmath.sqrt(y2_11)
    x2_11_sol2 = -x2_11_sol1
    print("The two possible values for the first coordinate of X2, (X2)_11, are:")
    print(f"  1. {x2_11_sol1}")
    print(f"  2. {x2_11_sol2}\n")

    # --- Final Sum ---
    # There are 2*2=4 solutions to the system (X1, X2).
    # We sum the first coordinate of the solution, which we interpret as (X1)_11.
    # The four solutions have (X1)_11 values: x1_11_sol1, x1_11_sol1, x1_11_sol2, x1_11_sol2.
    
    sum_of_coords = x1_11_sol1 + x1_11_sol1 + x1_11_sol2 + x1_11_sol2
    
    print("--- Calculating the Sum ---")
    print("There are 4 solutions to the system. The sum of their first coordinates, (X1)_11, is:")
    print("The final equation for the sum is:")
    print(f"({x1_11_sol1}) + ({x1_11_sol1}) + ({x1_11_sol2}) + ({x1_11_sol2}) = {sum_of_coords}")
    
    # The result is a complex number, but its value is 0. We display the real part.
    final_answer = sum_of_coords.real
    print(f"\nFinal Sum = {final_answer}")

solve_and_sum()