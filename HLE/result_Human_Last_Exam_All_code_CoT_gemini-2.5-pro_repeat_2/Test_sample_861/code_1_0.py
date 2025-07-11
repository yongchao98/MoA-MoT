from fractions import Fraction
import numpy as np

def solve_difference_equation():
    """
    Solves the difference equation 8y[n] - 6y[n-1] + y[n-2] = 1
    and calculates the value of E/A + (D*C)/B.
    """
    # Step 1 & 2: Find the homogeneous solution and roots B and D
    # The characteristic equation is 8r^2 - 6r + 1 = 0.
    # We solve for the roots r.
    char_coeffs = [8, -6, 1]
    roots = np.roots(char_coeffs)
    
    # By convention, assign the larger root to B and the smaller to D
    B_val = max(roots)
    D_val = min(roots)
    
    B = Fraction(B_val).limit_denominator()
    D = Fraction(D_val).limit_denominator()

    # Step 3: Find the particular solution E
    # For a constant input x[n]=1, we assume a constant particular solution y_p[n] = E.
    # 8E - 6E + E = 1 => 3E = 1
    E = Fraction(1, 3)

    # Step 4 & 5: Use initial conditions y[0]=1, y[-1]=2 to find A and C
    # The general solution is y[n] = A*(B**n) + C*(D**n) + E
    # For n=0: y[0] = A + C + E = 1  => A + C = 1 - E
    # For n=-1: y[-1] = A*(B**-1) + C*(D**-1) + E = 2 => A/B + C/D = 2 - E
    
    # System of linear equations:
    # eq1: A + C = 1 - E
    # eq2: A/B + C/D = 2 - E
    
    # We can represent this as M * [A, C]^T = V
    M = np.array([
        [1, 1],
        [1/B, 1/D]
    ])
    
    V = np.array([
        1 - E,
        2 - E
    ])
    
    # Solve for [A, C]
    coeffs = np.linalg.solve(M, V)
    A = Fraction(coeffs[0]).limit_denominator()
    C = Fraction(coeffs[1]).limit_denominator()

    # Step 6: Print the final equation and calculate the result
    def format_fraction(f):
        if f.denominator == 1:
            return str(f.numerator)
        return f"({f.numerator}/{f.denominator})"

    print(f"The closed-form solution is:")
    print(f"y[n] = {format_fraction(A)} * ({format_fraction(B)})^n + {format_fraction(C)} * ({format_fraction(D)})^n + {format_fraction(E)}\n")
    
    print("The values of the constants are:")
    print(f"A = {format_fraction(A)}")
    print(f"B = {format_fraction(B)}")
    print(f"C = {format_fraction(C)}")
    print(f"D = {format_fraction(D)}")
    print(f"E = {format_fraction(E)}\n")
    
    # Calculate the final expression: E/A + (D*C)/B
    result = E/A + (D*C)/B
    
    print(f"The value of E/A + (D*C)/B is {format_fraction(result)}, which is {float(result)}.")
    
    return float(result)

# Execute the function to get the final answer
final_answer = solve_difference_equation()
print(f'<<<{final_answer}>>>')
