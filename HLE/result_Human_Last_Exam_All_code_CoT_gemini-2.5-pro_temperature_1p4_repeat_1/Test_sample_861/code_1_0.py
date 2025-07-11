import numpy as np
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the required expression.
    """
    # Step 1: Find the homogeneous solution roots
    # The characteristic equation is 8r^2 - 6r + 1 = 0
    char_coeffs = [8, -6, 1]
    roots = np.roots(char_coeffs)
    
    # By convention, we assign the larger root to B and the smaller to D.
    roots = sorted(roots, reverse=True)
    B_val, D_val = roots[0], roots[1]

    # Use Fraction for precision
    B = Fraction(B_val).limit_denominator()
    D = Fraction(D_val).limit_denominator()

    # Step 2: Find the particular solution
    # We assume a particular solution y_p[n] = E (a constant).
    # Substituting into the original equation: 8E - 6E + E = 1 => 3E = 1
    E = Fraction(1, 3)

    # Step 3: Use initial conditions y[0]=1 and y[-1]=2 to find A and C
    # The general solution is y[n] = A*B^n + C*D^n + E
    # For n=0:  y[0] = A*B^0 + C*D^0 + E = A + C + E = 1
    # For n=-1: y[-1] = A*B^-1 + C*D^-1 + E = A/B + C/D + E = 2

    # This gives a system of two linear equations for A and C:
    # 1) A + C = 1 - E
    # 2) A/B + C/D = 2 - E
    
    # We can write this in matrix form: M * [A, C]^T = V^T
    M = np.array([
        [1, 1],
        [1/B, 1/D]
    ], dtype=float)
    
    V = np.array([
        1 - E,
        2 - E
    ], dtype=float)
    
    # Step 4: Solve for A and C
    try:
        solution = np.linalg.solve(M, V)
        A_val, C_val = solution[0], solution[1]
    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")
        return

    # Convert results back to fractions for precise representation
    A = Fraction(A_val).limit_denominator()
    C = Fraction(C_val).limit_denominator()

    # Step 5: Display the full solution
    print("The closed-form solution is of the form: y[n] = A*(B)^n + C*(D)^n + E")
    print("\nThe determined values for the constants are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    
    print("\nThus, the specific solution is:")
    print(f"y[n] = {A} * ({B})^n + {C} * ({D})^n + {E}")

    # Step 6: Calculate the final expression
    result_expr = E/A + (D*C)/B
    
    print("\nFinally, we calculate the expression E/A + (D*C)/B:")
    print(f"E/A = ({E}) / ({A}) = {E/A}")
    print(f"D*C = ({D}) * ({C}) = {D*C}")
    print(f"(D*C)/B = ({D*C}) / ({B}) = {(D*C)/B}")
    print(f"Result = E/A + (D*C)/B = {E/A} + {(D*C)/B} = {result_expr}")
    
    # The final numerical value
    print(f"\nThe final value is {float(result_expr)}")

solve_difference_equation()