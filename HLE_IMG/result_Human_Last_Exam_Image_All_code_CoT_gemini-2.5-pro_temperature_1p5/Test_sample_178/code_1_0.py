import numpy as np

def solve_temperature_distribution():
    """
    Calculates the steady-state temperature distribution on a square plate
    by solving the Poisson equation using the finite difference method.
    """
    # --- 1. Problem Definition & Grid Parameters ---
    num_intervals = 3
    plate_size = 0.3
    dx = plate_size / num_intervals
    dy = plate_size / num_intervals

    # Mapping of unknown labels to their (i, j) indices
    # (i corresponds to x, j corresponds to y)
    unknowns_map = {
        'T1': (1, 2), 'T2': (2, 2),
        'T3': (1, 1), 'T4': (2, 1)
    }
    unknowns_order = ['T1', 'T2', 'T3', 'T4']

    # --- 2. Calculate Coefficients ---
    lambda_val = (dx / dy)**2
    alpha_val = lambda_val + 1

    print("Step 1: Calculate coefficients")
    print(f"The grid has {num_intervals} intervals, so Δx = {dx:.1f} and Δy = {dy:.1f}.")
    print(f"λ = (Δx/Δy)² = ({dx:.1f}/{dy:.1f})² = {lambda_val:.2f}")
    print(f"α = λ + 1 = {lambda_val:.2f} + 1 = {alpha_val:.2f}\n")

    # --- 3. Set up the System of Equations AT = B ---
    print("Step 2: Formulate the system of linear equations.")
    # The finite difference equation is:
    # 2*α*T(i,j) - λ*[T(i,j+1)+T(i,j-1)] - [T(i+1,j)+T(i-1,j)] = -dx²*f(x,y)
    # For α=2, λ=1, it becomes: 4*T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -x*y

    A = np.zeros((4, 4))
    B = np.zeros(4)

    # Boundary Conditions
    T_left = 0.0      # T(0, y) = 0
    T_right = 1.0     # T(0.3, y) = 1
    T_bottom = 0.0    # T(x, 0) = 0
    T_top = 0.5       # T(x, 0.3) = 0.5

    # Equation for T1 (i=1, j=2) -> x=0.1, y=0.2
    # 4*T1 - T(top) - T(bottom_neighbor=T3) - T(right_neighbor=T2) - T(left) = -x*y
    A[0, 0] = 4  # T1
    A[0, 1] = -1 # T2
    A[0, 2] = -1 # T3
    B[0] = T_top + T_left - (0.1 * 0.2)

    # Equation for T2 (i=2, j=2) -> x=0.2, y=0.2
    # 4*T2 - T(top) - T(bottom_neighbor=T4) - T(right) - T(left_neighbor=T1) = -x*y
    A[1, 1] = 4  # T2
    A[1, 0] = -1 # T1
    A[1, 3] = -1 # T4
    B[1] = T_top + T_right - (0.2 * 0.2)

    # Equation for T3 (i=1, j=1) -> x=0.1, y=0.1
    # 4*T3 - T(top_neighbor=T1) - T(bottom) - T(right_neighbor=T4) - T(left) = -x*y
    A[2, 2] = 4  # T3
    A[2, 0] = -1 # T1
    A[2, 3] = -1 # T4
    B[2] = T_bottom + T_left - (0.1 * 0.1)

    # Equation for T4 (i=2, j=1) -> x=0.2, y=0.1
    # 4*T4 - T(top_neighbor=T2) - T(bottom) - T(right) - T(left_neighbor=T3) = -x*y
    A[3, 3] = 4  # T4
    A[3, 1] = -1 # T2
    A[3, 2] = -1 # T3
    B[3] = T_bottom + T_right - (0.2 * 0.1)
    
    print("The equations are formed by isolating unknown temperatures on the left side:")
    for i in range(4):
        eq_parts = []
        for j in range(4):
            if A[i, j] != 0:
                eq_parts.append(f"{A[i, j]:+.2f} * {unknowns_order[j]}")
        print(f"  { ' '.join(eq_parts) } = {B[i]:.4f}")
    
    # --- 4. Solve the System ---
    print("\nStep 3: Solve the system AT = B.")
    try:
        T_vector = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("Error: The matrix A is singular, the system cannot be solved.")
        return

    # --- 5. Output the Result ---
    print("\nThe vector value for {T} is:")
    print(f"{{T}} = [T1, T2, T3, T4] = [{T_vector[0]:.4f}, {T_vector[1]:.4f}, {T_vector[2]:.4f}, {T_vector[3]:.4f}]")

solve_temperature_distribution()