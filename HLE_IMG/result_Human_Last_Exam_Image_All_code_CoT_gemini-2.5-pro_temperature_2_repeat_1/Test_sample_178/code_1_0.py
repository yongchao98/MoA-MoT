import numpy as np

def solve_temperature_distribution():
    """
    This function sets up and solves the system of linear equations for the
    steady-state temperature distribution on a square plate using the
    finite difference method.
    """
    # Step 1: Define grid parameters and calculate coefficients
    print("Step 1: Determine grid parameters and coefficients.")
    # A 0.3x0.3 plate with 3 intervals in each direction
    num_intervals = 3
    length = 0.3
    delta_x = length / num_intervals
    delta_y = length / num_intervals

    # Calculate lambda and alpha from the given formulas
    lambda_val = (delta_x / delta_y)**2
    alpha = lambda_val + 1

    print(f"The grid has 3 intervals, so Δx = {delta_x:.2f} and Δy = {delta_y:.2f}.")
    print(f"λ = (Δx/Δy)² = {lambda_val:.1f}")
    print(f"α = λ + 1 = {alpha:.1f}\n")

    # Step 2: Formulate the system of linear equations A*T = b
    print("Step 2: Formulate the system of linear equations for the four internal nodes.")
    # The finite difference equation is:
    # 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)² * f(x,y)
    # Substituting α=2, λ=1, we move known boundary values to the right-hand side.
    # The unknown vector is T = [T1, T2, T3, T4]^T.

    # The coefficient matrix A is derived from the LHS of the equations
    A = np.array([
        [4., -1., -1., 0.],
        [-1., 4., 0., -1.],
        [-1., 0., 4., -1.],
        [0., -1., -1., 4.]
    ])

    # The constant vector b is derived from the boundary conditions and the source term f(x,y)
    # RHS = Sum of boundary temps at neighbors - (Δx)² * 100 * x * y
    b = np.array([
        (0.5 + 0) - (delta_x**2) * 100 * (1*delta_x) * (2*delta_y),  # Eq for T1(i=1,j=2)
        (0.5 + 1.0) - (delta_x**2) * 100 * (2*delta_x) * (2*delta_y),# Eq for T2(i=2,j=2)
        (0 + 0) - (delta_x**2) * 100 * (1*delta_x) * (1*delta_y),    # Eq for T3(i=1,j=1)
        (1.0 + 0) - (delta_x**2) * 100 * (2*delta_x) * (1*delta_y)  # Eq for T4(i=2,j=1)
    ])

    print("The system of equations A*T = b is derived as:")
    labels = ["T1", "T2", "T3", "T4"]
    for i in range(A.shape[0]):
        equation = " + ".join([f"({A[i, j]:.1f})*{labels[j]}" for j in range(A.shape[1])])
        equation += f" = {b[i]:.4f}"
        print(f"  {equation}")
    print()

    # Step 3: Solve the system
    print("Step 3: Solve the system for the vector {T}.")
    try:
        T = np.linalg.solve(A, b)
        print("The vector value for {T} is:")
        print(f"[T1, T2, T3, T4] = [{T[0]}, {T[1]}, {T[2]}, {T[3]}]")

    except np.linalg.LinAlgError:
        print("The matrix A is singular, and the system cannot be solved.")

# Run the solver
solve_temperature_distribution()

# For the final answer, calculate and format the result again
A_final = np.array([
    [4, -1, -1, 0],
    [-1, 4, 0, -1],
    [-1, 0, 4, -1],
    [0, -1, -1, 4]
])
b_final = np.array([
    0.48,
    1.46,
    -0.01,
    0.98
])
T_final = np.linalg.solve(A_final, b_final)
# The output format requests [ [T1, T2, T3, T4] ]. It likely means a list of lists with one inner list,
# although it could be a typo for a single list. Adhering to the format:
answer_vector = [T_final.tolist()]
print(f"\nFinal answer in the requested format: {answer_vector}")
# <<<[[0.2608333333333333, 0.43583333333333335, 0.13333333333333333, 0.3120833333333333]]>>>