import numpy as np

def solve_bvp():
    """
    Solves the boundary value problem and calculates the final expression.
    """
    # 1. Define constants
    sqrt5 = np.sqrt(5)
    lambda1 = (1 + sqrt5) / 2
    lambda2 = (1 - sqrt5) / 2

    # 2. Define the boundary points
    t0 = 0
    t1 = np.log(5)

    # 3. Set up the linear system for the coefficients A and B of d(t)
    # d(t) = A*exp(lambda1*t) + B*exp(lambda2*t) - 2/sqrt(5)
    # d(t0) = 0  => A + B = 2/sqrt(5)
    # d(t1) = 0  => A*exp(lambda1*t1) + B*exp(lambda2*t1) = 2/sqrt(5)
    
    # Matrix for the system [c0, c1] for [A, B]
    # c0*A + c1*B = rhs
    # [ 1,                  1                ] [A] = [2/sqrt5]
    # [ exp(lambda1*t1),    exp(lambda2*t1)    ] [B] = [2/sqrt5]
    
    # Use numpy to solve for A and B
    a_matrix = np.array([
        [1, 1],
        [np.exp(lambda1 * t1), np.exp(lambda2 * t1)]
    ])
    b_vector = np.array([2 / sqrt5, 2 / sqrt5])
    
    try:
        coeffs_A_B = np.linalg.solve(a_matrix, b_vector)
        A = coeffs_A_B[0]
        B = coeffs_A_B[1]
    except np.linalg.LinAlgError:
        print("Matrix is singular. Could not solve for coefficients A and B.")
        return

    # 4. Define the solution phi_0(t)
    # phi_0(t) = (A - 1/sqrt5)*exp(lambda1*t) + (B - 1/sqrt5)*exp(lambda2*t)
    C1 = A - 1/sqrt5
    C2 = B - 1/sqrt5
    
    def phi_0(t):
        return C1 * np.exp(lambda1 * t) + C2 * np.exp(lambda2 * t)

    # 5. Calculate the final value
    t_f = np.log(10**10)
    
    # The value is -phi_0(t_f) + 2/sqrt(5)
    value = -phi_0(t_f) + 2/sqrt5
    
    # Output the result step-by-step
    print(f"The golden ratio roots are:")
    print(f"λ₁ = (1 + √5)/2 = {lambda1}")
    print(f"λ₂ = (1 - √5)/2 = {lambda2}\n")
    
    print(f"The solution for φ₀(t) has the form C₁*exp(λ₁*t) + C₂*exp(λ₂*t).")
    print(f"Solving the boundary conditions gives the constants:")
    print(f"C₁ = {C1}")
    print(f"C₂ = {C2}\n")

    phi_val_tf = phi_0(t_f)
    print(f"The value to evaluate is t_f = ln(10¹⁰) = {t_f}")
    print(f"φ₀(t_f) = {phi_val_tf}\n")
    
    final_expr_val_1 = -phi_val_tf
    final_expr_val_2 = 2/sqrt5
    print(f"The expression to calculate is: -φ₀(t_f) + 2/√5")
    print(f"Which is: -({phi_val_tf}) + {final_expr_val_2}")
    print(f"Result = {value}")

solve_bvp()