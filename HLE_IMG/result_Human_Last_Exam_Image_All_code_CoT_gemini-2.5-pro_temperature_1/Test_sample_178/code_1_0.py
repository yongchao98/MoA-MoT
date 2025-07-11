import numpy as np

def solve_temperature_distribution():
    """
    Calculates the steady-state temperature distribution on a square plate
    by solving the Poisson equation using the finite difference method.
    """
    # --- Step 1: Define problem parameters and calculate coefficients ---
    # The plate is 0.3x0.3 with 3 intervals, so step sizes are 0.1.
    delta_x = 0.1
    delta_y = 0.1

    # Calculate lambda and alpha as per the given formulas
    lambda_val = (delta_x / delta_y)**2
    alpha = lambda_val + 1

    print("Step 1: Calculated Coefficients")
    print(f"λ = (Δx/Δy)² = ({delta_x}/{delta_y})² = {lambda_val:.2f}")
    print(f"α = λ + 1 = {lambda_val:.2f} + 1 = {alpha:.2f}\n")

    # --- Step 2: Set up the system of linear equations (A * T = B) ---
    # The finite difference equation is:
    # 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)² * f(x(i), y(i))
    # With α=2, λ=1, f(x,y)=100xy, and Δx=0.1, the equation simplifies to:
    # 4*T(i,j) - T_top - T_bottom - T_right - T_left = -x*y
    
    # Let the unknown vector be {T} = [T1, T2, T3, T4]^T.
    # The equations are derived by applying the formula at each of the 4 interior nodes.
    # T1 at (x=0.1, y=0.2), T2 at (x=0.2, y=0.2)
    # T3 at (x=0.1, y=0.1), T4 at (x=0.2, y=0.1)

    # Boundary conditions:
    # T(0,y)=0, T(0.3,y)=1, T(x,0)=0, T(x,0.3)=0.5

    # Equation for T1(x=0.1, y=0.2): 4*T1 - T2 - T3 = 0.5 - (0.1*0.2) = 0.48
    # Equation for T2(x=0.2, y=0.2): 4*T2 - T1 - T4 = 1 + 0.5 - (0.2*0.2) = 1.46
    # Equation for T3(x=0.1, y=0.1): 4*T3 - T1 - T4 = 0 + 0 - (0.1*0.1) = -0.01
    # Equation for T4(x=0.2, y=0.1): 4*T4 - T2 - T3 = 1 + 0 - (0.2*0.1) = 0.98

    # Coefficient matrix A
    A = np.array([
        [ 4., -1., -1.,  0.],
        [-1.,  4.,  0., -1.],
        [-1.,  0.,  4., -1.],
        [ 0., -1., -1.,  4.]
    ])

    # Constant vector B
    B = np.array([0.48, 1.46, -0.01, 0.98])

    # --- Step 3: Print the system of equations ---
    print("Step 2 & 3: The System of Linear Equations (A * {T} = {B})")
    print("The final equations for T1, T2, T3, and T4 are:")
    
    # We will print each equation in a readable format.
    for i in range(A.shape[0]):
        equation_str = []
        # First term
        if A[i, 0] != 0:
            equation_str.append(f"{A[i, 0]:.2f}*T1")
        # Subsequent terms
        for j in range(1, A.shape[1]):
            if A[i, j] != 0:
                sign = " + " if A[i, j] > 0 else " - "
                equation_str.append(f"{sign}{abs(A[i, j]):.2f}*T{j+1}")
        
        print(f"{''.join(equation_str)} = {B[i]:.2f}")

    # --- Step 4: Solve the system for {T} ---
    try:
        T = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        print("\nError: The system of equations could not be solved.")
        return

    # --- Step 5: Display the final result ---
    print("\nStep 4 & 5: Solving for {T} and Final Answer")
    print("The vector value for {T} is:")
    
    # Format the output to be clean
    formatted_T = [f"{val:.4f}" for val in T]
    print(f"{{T}} = [{', '.join(formatted_T)}]")

solve_temperature_distribution()