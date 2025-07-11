import numpy as np
from fractions import Fraction

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Step 1: Find the homogeneous solution by solving the characteristic equation
    # 8r^2 - 6r + 1 = 0
    char_eq_coeffs = [8, -6, 1]
    roots = np.roots(char_eq_coeffs)

    # By convention, assign the larger root to B and the smaller to D
    B = Fraction(max(roots)).limit_denominator()
    D = Fraction(min(roots)).limit_denominator()

    # Step 2: Find the particular solution
    # For a constant input of 1, assume y_p[n] = K
    # 8K - 6K + K = 1  => 3K = 1
    E = Fraction(1, 3)

    # Step 3 & 4: Use initial conditions to find A and C
    # The general solution is y[n] = A*B**n + C*D**n + E
    # Using y[0] = 1: A*B**0 + C*D**0 + E = 1 => A + C = 1 - E
    # Using y[-1] = 2: A*B**-1 + C*D**-1 + E = 2 => A/B + C/D = 2 - E

    # Setup the system of linear equations to solve for A and C:
    # Eq1: 1*A + 1*C = 1 - E
    # Eq2: (1/B)*A + (1/D)*C = 2 - E
    matrix_coeffs = np.array([
        [1, 1],
        [1/B, 1/D]
    ], dtype=float)

    constants = np.array([
        1 - E,
        2 - E
    ], dtype=float)

    # Solve for [A, C]
    solution = np.linalg.solve(matrix_coeffs, constants)
    A = Fraction(solution[0]).limit_denominator()
    C = Fraction(solution[1]).limit_denominator()

    # Step 5: Calculate the final expression E/A + (D*C)/B
    result = E/A + (D*C)/B

    # Output the results
    print(f"The difference equation is: 8y[n] - 6y[n-1] + y[n-2] = 1")
    print(f"The general form of the solution is: y[n] = A * B^n + C * D^n + E")
    print("\nSolving for the coefficients:")
    print(f"  The roots of the characteristic equation 8r^2 - 6r + 1 = 0 are {max(roots):.2f} and {min(roots):.2f}.")
    print(f"  Assigning the larger root to B and the smaller to D:")
    print(f"  B = {B}")
    print(f"  D = {D}")
    print(f"  The particular solution gives E:")
    print(f"  E = {E}")
    print(f"  Solving for A and C using initial conditions y[0]=1, y[-1]=2:")
    print(f"  A = {A}")
    print(f"  C = {C}")

    print("\nThe final solution is:")
    print(f"y[n] = ({A}) * ({B})^n + ({C}) * ({D})^n + ({E})")

    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"({E}) / ({A}) + (({D}) * ({C})) / ({B})")
    # Show intermediate steps of the calculation
    term1 = E/A
    term2_numerator = D*C
    term2 = term2_numerator / B
    print(f"= ({term1}) + ({term2_numerator}) / ({B})")
    print(f"= ({term1}) + ({term2})")
    print(f"= {result}")

    print(f"\nThe final answer is {float(result)}")

if __name__ == '__main__':
    solve_difference_equation()