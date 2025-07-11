import numpy as np

def solve_temperature_distribution():
    """
    Calculates the steady-state temperature distribution on a square plate
    by solving the Poisson equation using the finite difference method.
    """
    # Step 1: Define grid parameters and constants
    num_intervals = 3
    length = 0.3
    delta_x = length / num_intervals
    delta_y = length / num_intervals

    lambda_val = (delta_x / delta_y)**2
    alpha = lambda_val + 1

    print(f"Calculated Parameters:")
    print(f"Δx = {delta_x:.1f}")
    print(f"Δy = {delta_y:.1f}")
    print(f"λ = (Δx/Δy)^2 = {lambda_val:.1f}")
    print(f"α = λ + 1 = {alpha:.1f}\n")

    # Step 2: Set up the system of linear equations A * T = B
    # The finite difference equation is:
    # 2*α*T - λ*(T_up + T_down) - (T_right + T_left) = - (Δx)^2 * f(x, y)
    # With α=2, λ=1, Δx=0.1, f(x,y)=100xy, it simplifies to:
    # 4*T - T_up - T_down - T_right - T_left = -x*y
    
    # Unknown vector is T = [T1, T2, T3, T4]^T
    # Boundary Conditions: T_left=0, T_right=1.0, T_bottom=0, T_top=0.5
    
    # Equation for T1 (x=0.1, y=0.2):
    # 4*T1 - T_top - T3 - T2 - T_left = -0.1*0.2
    # 4*T1 - 0.5 - T3 - T2 - 0 = -0.02 => 4*T1 - T2 - T3 = 0.48

    # Equation for T2 (x=0.2, y=0.2):
    # 4*T2 - T_top - T4 - T_right - T1 = -0.2*0.2
    # 4*T2 - 0.5 - T4 - 1.0 - T1 = -0.04 => -T1 + 4*T2 - T4 = 1.46

    # Equation for T3 (x=0.1, y=0.1):
    # 4*T3 - T1 - T_bottom - T4 - T_left = -0.1*0.1
    # 4*T3 - T1 - 0 - T4 - 0 = -0.01 => -T1 + 4*T3 - T4 = -0.01

    # Equation for T4 (x=0.2, y=0.1):
    # 4*T4 - T2 - T_bottom - T_right - T3 = -0.2*0.1
    # 4*T4 - T2 - 0 - 1.0 - T3 = -0.02 => -T2 - T3 + 4*T4 = 0.98

    # Coefficient matrix A
    A = np.array([
        [ 4, -1, -1,  0],
        [-1,  4,  0, -1],
        [-1,  0,  4, -1],
        [ 0, -1, -1,  4]
    ])

    # Right-hand side vector B
    B = np.array([0.48, 1.46, -0.01, 0.98])
    
    print("The system of linear equations is AT = B, where T = [T1, T2, T3, T4]^T.")
    print("The final equations are:")
    for i in range(len(B)):
        equation_str = []
        for j in range(len(B)):
            if A[i,j] != 0:
                equation_str.append(f"({A[i,j]:.1f})*T{j+1}")
        print(f"{' + '.join(equation_str).replace('+ (-', '- (')} = {B[i]:.2f}")
    
    # Step 3: Solve the system for T
    T = np.linalg.solve(A, B)

    # Step 4: Print the result in the required format
    print("\nThe vector value for {T} is:")
    print(f"[{T[0]:.4f}, {T[1]:.4f}, {T[2]:.4f}, {T[3]:.4f}]")

if __name__ == '__main__':
    solve_temperature_distribution()