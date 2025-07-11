import numpy as np

def solve_temperature_distribution():
    """
    Solves the steady-state temperature distribution on a square plate
    using the finite difference method for the Poisson equation.
    """
    # Step 1: Discretize the domain and define parameters
    num_intervals = 3
    plate_size = 0.3
    delta_x = plate_size / num_intervals
    delta_y = plate_size / num_intervals

    # Boundary conditions
    T_right = 1.0  # T(0.3, y)
    T_top = 0.5    # T(x, 0.3)
    T_left = 0.0   # T(0, y)
    T_bottom = 0.0 # T(x, 0)

    # Forcing function f(x,y) = 100xy
    def f(x, y):
        return 100 * x * y

    # Step 2: Calculate coefficients lambda and alpha
    lam = (delta_x / delta_y)**2
    alpha = lam + 1

    print("Step 1: Calculate coefficients.")
    print(f"The grid spacing is delta_x = {delta_x:.1f} and delta_y = {delta_y:.1f}.")
    print(f"λ = (Δx/Δy)² = ({delta_x:.1f}/{delta_y:.1f})² = {lam:.1f}")
    print(f"α = λ + 1 = {lam:.1f} + 1 = {alpha:.1f}")
    print("-" * 40)

    # Step 3: Formulate the system of linear equations A*T = b
    # The general equation is: 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)²*f(x,y)
    # Substituting α=2, λ=1: 4*T(i,j) - [T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)²*f(x,y)
    
    # Coordinates of unknown points:
    # T1 at (x=0.1, y=0.2), T2 at (x=0.2, y=0.2)
    # T3 at (x=0.1, y=0.1), T4 at (x=0.2, y=0.1)

    # Calculate the right-hand side vector b, moving known boundary terms to the right
    # b_node = sum(known_neighbor_T) - (delta_x**2) * f(x,y)
    b1 = T_top + T_left - (delta_x**2) * f(0.1, 0.2)
    b2 = T_top + T_right - (delta_x**2) * f(0.2, 0.2)
    b3 = T_bottom + T_left - (delta_x**2) * f(0.1, 0.1)
    b4 = T_bottom + T_right - (delta_x**2) * f(0.2, 0.1)
    
    # Define matrix A based on connectivity of unknown nodes
    # A's coefficients are -1 for neighboring unknowns and 4 on the diagonal.
    # T vector is [T1, T2, T3, T4]
    A = np.array([
        [4, -1, -1,  0],  # Eq for T1: Neighbors are T_top, T_left, T2, T3
        [-1, 4,  0, -1],  # Eq for T2: Neighbors are T_top, T_right, T1, T4
        [-1, 0,  4, -1],  # Eq for T3: Neighbors are T1, T4, T_left, T_bottom
        [0, -1, -1,  4]   # Eq for T4: Neighbors are T2, T3, T_right, T_bottom
    ])

    b = np.array([b1, b2, b3, b4])

    print("Step 2: Formulate the system of equations A*T = b.")
    print("The system is derived by applying the finite difference formula to each unknown node.")
    
    print("\nEquation for T1:")
    print(f"4*T1 - T2 - T3 - T_top - T_left = -({delta_x})² * f(0.1, 0.2)")
    print(f"4*T1 - T2 - T3 - {T_top} - {T_left} = -({delta_x**2:.2f}) * ({f(0.1, 0.2):.1f})")
    print(f"Final form: ({A[0,0]})*T1 + ({A[0,1]})*T2 + ({A[0,2]})*T3 + ({A[0,3]})*T4 = {b[0]:.2f}")

    print("\nEquation for T2:")
    print(f"4*T2 - T1 - T4 - T_top - T_right = -({delta_x})² * f(0.2, 0.2)")
    print(f"4*T2 - T1 - T4 - {T_top} - {T_right} = -({delta_x**2:.2f}) * ({f(0.2, 0.2):.1f})")
    print(f"Final form: ({A[1,0]})*T1 + ({A[1,1]})*T2 + ({A[1,2]})*T3 + ({A[1,3]})*T4 = {b[1]:.2f}")

    print("\nEquation for T3:")
    print(f"4*T3 - T1 - T4 - T_bottom - T_left = -({delta_x})² * f(0.1, 0.1)")
    print(f"4*T3 - T1 - T4 - {T_bottom} - {T_left} = -({delta_x**2:.2f}) * ({f(0.1, 0.1):.1f})")
    print(f"Final form: ({A[2,0]})*T1 + ({A[2,1]})*T2 + ({A[2,2]})*T3 + ({A[2,3]})*T4 = {b[2]:.2f}")

    print("\nEquation for T4:")
    print(f"4*T4 - T2 - T3 - T_bottom - T_right = -({delta_x})² * f(0.2, 0.1)")
    print(f"4*T4 - T2 - T3 - {T_bottom} - {T_right} = -({delta_x**2:.2f}) * ({f(0.2, 0.1):.1f})")
    print(f"Final form: ({A[3,0]})*T1 + ({A[3,1]})*T2 + ({A[3,2]})*T3 + ({A[3,3]})*T4 = {b[3]:.2f}")
    print("-" * 40)
    
    # Step 4: Solve the system
    T = np.linalg.solve(A, b)
    
    # Step 5: Display the result
    print("Step 3: Solve for {T} = [T1, T2, T3, T4].")
    print("\nThe vector value for {T} is:")
    print(f"{T.tolist()}")

solve_temperature_distribution()
<<<[0.2642857142857143, 0.4914285714285714, 0.17428571428571428, 0.4114285714285714]>>>