import numpy as np

def solve_difference_equation():
    """
    Solves the given difference equation and calculates the final expression.
    """
    # Step 1 & 2: Find roots of characteristic equation and the particular solution
    # Homogeneous equation: 8y[n] - 6y[n-1] + y[n-2] = 0
    # Characteristic equation: 8r^2 - 6r + 1 = 0
    char_eq_coeffs = [8, -6, 1]
    roots = np.roots(char_eq_coeffs)

    # Particular solution: 8K - 6K + K = 1  => 3K = 1
    E = 1/3

    # Step 3 & 4: Use initial conditions to find coefficients
    # y[0] = 1, y[-1] = 2
    # Convention: B is the larger root, D is the smaller root.
    B = max(roots)
    D = min(roots)

    # General solution form: y[n] = A*(B)^n + C*(D)^n + E
    # At n=0: A + C + E = y[0] => A + C = y[0] - E
    # At n=-1: A*B^-1 + C*D^-1 + E = y[-1] => A/B + C/D = y[-1] - E
    y0 = 1
    ym1 = 2
    
    # System of linear equations:
    # 1*A + 1*C = y0 - E
    # (1/B)*A + (1/D)*C = ym1 - E
    matrix_coeffs = np.array([[1, 1], [1/B, 1/D]])
    constants_vector = np.array([y0 - E, ym1 - E])
    
    try:
        solution = np.linalg.solve(matrix_coeffs, constants_vector)
        A = solution[0]
        C = solution[1]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # Step 5 & 6: Print results and calculate final expression
    print("The closed-form solution is y[n] = A*(B)^n + C*(D)^n + E")
    print(f"Based on the convention B > D:")
    print(f"A = {A:.4f}")
    print(f"B = {B:.4f}")
    print(f"C = {C:.4f}")
    print(f"D = {D:.4f}")
    print(f"E = {E:.4f}")
    print("\nThe full equation is:")
    print(f"y[n] = ({A:.4f}) * ({B:.4f})^n + ({C:.4f}) * ({D:.4f})^n + {E:.4f}")

    # Calculate the final value
    result = E/A + (D*C)/B
    
    print("\nCalculating the expression E/A + (D*C)/B:")
    print(f"E/A = ({E:.4f}) / ({A:.4f})")
    print(f"(D*C)/B = (({D:.4f}) * ({C:.4f})) / ({B:.4f})")
    print(f"Result = {result:.4f}")

    # Return the precise value for the final answer tag
    return result

if __name__ == '__main__':
    final_value = solve_difference_equation()
    # The final answer is requested in a specific format at the end
    # We will print the final numerical value wrapped in <<<>>>
    # Let's use a higher precision for the final tag
    print(f"\nFinal answer rounded to 4 decimal places: {final_value:.4f}")

# The user wants the value of E/A + (D*C)/B
# Using fractions to get the exact value:
A_f, B_f, C_f, D_f, E_f = (1/2, 1/2, 1/6, 1/4, 1/3)
final_exact_value = (E_f / A_f) + (D_f * C_f) / B_f
# (1/3)/(1/2) + (1/4 * 1/6)/(1/2) = 2/3 + (1/24)/(1/2) = 2/3 + 1/12 = 8/12 + 1/12 = 9/12 = 0.75
print(f"<<<{final_exact_value}>>>")