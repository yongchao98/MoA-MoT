import sympy

def solve_poisson_equation():
    """
    Solves the finite difference approximation of the Poisson equation for the given grid and boundary conditions.
    """
    # Step 1: Define Grid Parameters
    num_intervals = 3
    length = 0.3
    dx = length / num_intervals
    dy = length / num_intervals

    # Step 2: Calculate Coefficients
    lambda_val = (dx / dy)**2
    alpha_val = lambda_val + 1

    print(f"Calculated Parameters:")
    print(f"Δx = {dx}")
    print(f"Δy = {dy}")
    print(f"λ = {lambda_val}")
    print(f"α = {alpha_val}\n")

    # Define coordinates for T1, T2, T3, T4
    coords = {
        'T1': (dx * 1, dy * 2),  # (0.1, 0.2)
        'T2': (dx * 2, dy * 2),  # (0.2, 0.2)
        'T3': (dx * 1, dy * 1),  # (0.1, 0.1)
        'T4': (dx * 2, dy * 1),  # (0.2, 0.1)
    }

    # Define boundary conditions
    T_x0 = 0      # T(0, y) = 0
    T_xL = 1      # T(0.3, y) = 1
    T_y0 = 0      # T(x, 0) = 0
    T_yL = 0.5    # T(x, 0.3) = 0.5
    
    # Define the source function f(x,y)
    def f(x, y):
        return 100 * x * y

    # Step 3 & 4: Formulate Linear Equations and Apply Boundary Conditions
    # The general equation is: 2*α*T(i,j) - λ*[T(i,j+1) + T(i,j-1)] - [T(i+1,j) + T(i-1,j)] = -(Δx)^2*f(x,y)
    # With λ=1, α=2, it simplifies to: 4T(i,j) - T(i,j+1) - T(i,j-1) - T(i+1,j) - T(i-1,j) = -dx^2*f(x,y)
    # We rearrange to put known boundary values on the right-hand side.
    
    # Equation for T1 at (0.1, 0.2)
    # Neighbors: T_yL (top), T3 (bottom), T2 (right), T_x0 (left)
    # 4*T1 - T_yL - T3 - T2 - T_x0 = -dx^2 * f(0.1, 0.2)
    # 4*T1 - T2 - T3 = T_yL + T_x0 - dx^2 * f(0.1, 0.2)
    b1 = T_yL + T_x0 - dx**2 * f(coords['T1'][0], coords['T1'][1])

    # Equation for T2 at (0.2, 0.2)
    # Neighbors: T_yL (top), T4 (bottom), T_xL (right), T1 (left)
    # 4*T2 - T_yL - T4 - T_xL - T1 = -dx^2 * f(0.2, 0.2)
    # -T1 + 4*T2 - T4 = T_yL + T_xL - dx^2 * f(0.2, 0.2)
    b2 = T_yL + T_xL - dx**2 * f(coords['T2'][0], coords['T2'][1])

    # Equation for T3 at (0.1, 0.1)
    # Neighbors: T1 (top), T_y0 (bottom), T4 (right), T_x0 (left)
    # 4*T3 - T1 - T_y0 - T4 - T_x0 = -dx^2 * f(0.1, 0.1)
    # -T1 + 4*T3 - T4 = T_y0 + T_x0 - dx^2 * f(0.1, 0.1)
    b3 = T_y0 + T_x0 - dx**2 * f(coords['T3'][0], coords['T3'][1])

    # Equation for T4 at (0.2, 0.1)
    # Neighbors: T2 (top), T_y0 (bottom), T_xL (right), T3 (left)
    # 4*T4 - T2 - T_y0 - T_xL - T3 = -dx^2 * f(0.2, 0.1)
    # -T2 - T3 + 4*T4 = T_y0 + T_xL - dx^2 * f(0.2, 0.1)
    b4 = T_y0 + T_xL - dx**2 * f(coords['T4'][0], coords['T4'][1])

    # Step 5: Construct the Matrix System A * {T} = {b}
    # Using exact fractions with sympy to avoid floating point errors
    A = sympy.Matrix([
        [4, -1, -1,  0],
        [-1, 4,  0, -1],
        [-1, 0,  4, -1],
        [0, -1, -1,  4]
    ])
    
    b = sympy.Matrix([
        sympy.nsimplify(b1),
        sympy.nsimplify(b2),
        sympy.nsimplify(b3),
        sympy.nsimplify(b4)
    ])

    print("System of Linear Equations (A * {T} = b):")
    vars = ['T1', 'T2', 'T3', 'T4']
    for i in range(A.rows):
        eq_str = " + ".join([f"{A[i, j]:.2f} * {vars[j]}" for j in range(A.cols)])
        print(f"Eq {i+1}: {eq_str} = {float(b[i]):.2f}")
    print()
    
    # Step 6: Solve the System
    T_solution = A.solve(b)

    # Step 7: Display the Results
    T_vector = [float(T) for T in T_solution]
    
    print("Final vector value for {T}:")
    print(f"[T1, T2, T3, T4] = {T_vector}")
    
    return T_vector

if __name__ == '__main__':
    solution = solve_poisson_equation()
    print(f"\n<<<{solution}>>>")
