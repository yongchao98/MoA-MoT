import numpy as np
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given linear difference equation and calculates the final expression.
    """
    # Step 1: Find the Homogeneous Solution by solving the characteristic equation 8r^2 - 6r + 1 = 0
    char_eq_coeffs = [8, -6, 1]
    roots = np.roots(char_eq_coeffs)

    # By convention, assign the larger root to B and the smaller to D
    B = max(roots)
    D = min(roots)

    # Step 2: Find the Particular Solution E
    # For a constant input '1', assume a constant output y_p[n] = E
    # 8E - 6E + E = 1 => 3E = 1
    E = 1/3

    # Step 3 & 4: Use initial conditions y[0]=1 and y[-1]=2 to find A and C
    # y[0] = 1 => A*B^0 + C*D^0 + E = 1 => A + C = 1 - E
    # y[-1] = 2 => A*B^-1 + C*D^-1 + E = 2 => A/B + C/D = 2 - E
    
    # We have a system of two linear equations:
    # 1*A + 1*C = 1 - E
    # (1/B)*A + (1/D)*C = 2 - E
    
    matrix_coeffs = np.array([
        [1, 1],
        [1/B, 1/D]
    ])
    
    constants = np.array([
        1 - E,
        2 - E
    ])
    
    # Solve the system for [A, C]
    solution = np.linalg.solve(matrix_coeffs, constants)
    A = solution[0]
    C = solution[1]

    # Display the results as fractions for clarity
    A_f = Fraction(A).limit_denominator()
    B_f = Fraction(B).limit_denominator()
    C_f = Fraction(C).limit_denominator()
    D_f = Fraction(D).limit_denominator()
    E_f = Fraction(E).limit_denominator()
    
    print("Derived coefficients for the equation y[n] = A*(B)^n + C*(D)^n + E:")
    print(f"A = {A_f}")
    print(f"B = {B_f}")
    print(f"C = {C_f}")
    print(f"D = {D_f}")
    print(f"E = {E_f}")
    print("-" * 30)

    # Print the final equation with each number
    print("The final solution is:")
    print(f"y[n] = ({A_f}) * ({B_f})^n + ({C_f}) * ({D_f})^n + ({E_f})")
    print("-" * 30)

    # Step 5: Calculate the final expression E/A + (D*C)/B
    result = E/A + (D*C)/B
    result_f = Fraction(result).limit_denominator()

    print("Calculating the expression E/A + (D*C)/B:")
    print(f"E/A = ({E_f}) / ({A_f}) = {Fraction(E/A).limit_denominator()}")
    print(f"D*C = ({D_f}) * ({C_f}) = {Fraction(D*C).limit_denominator()}")
    print(f"(D*C)/B = ({Fraction(D*C).limit_denominator()}) / ({B_f}) = {Fraction((D*C)/B).limit_denominator()}")
    print(f"Final Result = {Fraction(E/A).limit_denominator()} + {Fraction((D*C)/B).limit_denominator()} = {result_f}")
    print("-" * 30)
    print(f"The final numerical value is: {result}")


solve_difference_equation()