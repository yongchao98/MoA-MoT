import numpy as np

def solve_and_sum_first_coordinates():
    """
    Solves the given system of matrix equations and finds the sum of the
    first coordinate of all solutions.
    """
    # --- Equation 1 ---
    # (5 0; 0 -5) X1^2 + X1^2 (6 0; 0 6) = (-53/12 0; 0 0)
    # This simplifies to (A1 + 6I) * X1^2 = C1
    
    A1 = np.array([[5., 0.], [0., -5.]])
    k1 = 6.
    C1 = np.array([[-53./12., 0.], [0., 0.]])

    # Calculate Y1 = X1^2
    M1 = A1 + k1 * np.identity(2)
    inv_M1 = np.linalg.inv(M1)
    Y1 = inv_M1 @ C1

    # The (1,1) element of Y1, denoted y1_11, is Y1[0, 0].
    # The first coordinate of X1, let's call it x1_11, satisfies x1_11^2 = y1_11.
    # The solutions are sqrt(y1_11) and -sqrt(y1_11).
    y1_11 = Y1[0, 0]
    # Add 0j to ensure the argument is treated as a complex number for sqrt
    sol1_coord1 = np.sqrt(y1_11 + 0j)
    sol1_coord2 = -sol1_coord1

    # --- Equation 2 ---
    # (4 0; 0 -5) X2^2 + X2^2 (6 0; 0 6) = (-3/11 0; 0 0)
    # This simplifies to (A2 + 6I) * X2^2 = C2
    
    A2 = np.array([[4., 0.], [0., -5.]])
    k2 = 6.
    C2 = np.array([[-3./11., 0.], [0., 0.]])
    
    # Calculate Y2 = X2^2
    M2 = A2 + k2 * np.identity(2)
    inv_M2 = np.linalg.inv(M2)
    Y2 = inv_M2 @ C2

    # The first coordinate of X2, x2_11, satisfies x2_11^2 = Y2[0, 0].
    y2_11 = Y2[0, 0]
    sol2_coord1 = np.sqrt(y2_11 + 0j)
    sol2_coord2 = -sol2_coord1

    # The set of solutions includes all found matrices. The sum is over all their first coordinates.
    total_sum = sol1_coord1 + sol1_coord2 + sol2_coord1 + sol2_coord2
    
    # Print the calculation process and the final result as requested.
    print("The first coordinates of the solutions for X1 are:")
    print(f"Solution 1: {sol1_coord1}")
    print(f"Solution 2: {sol1_coord2}\n")
    
    print("The first coordinates of the solutions for X2 are:")
    print(f"Solution 3: {sol2_coord1}")
    print(f"Solution 4: {sol2_coord2}\n")
    
    print("The final equation is the sum of these coordinates:")
    print(f"Sum = ({sol1_coord1}) + ({sol1_coord2}) + ({sol2_coord1}) + ({sol2_coord2})")
    
    # The imaginary part will be zero due to cancellation.
    print(f"Final Sum = {total_sum.real}")

solve_and_sum_first_coordinates()