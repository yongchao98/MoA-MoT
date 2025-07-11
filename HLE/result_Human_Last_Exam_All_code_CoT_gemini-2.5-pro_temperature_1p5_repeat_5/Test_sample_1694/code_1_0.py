import numpy as np
import cmath

def solve_and_sum():
    """
    Solves the two matrix equations and finds the sum of the first coordinates of all solutions.
    """
    # --- First Equation ---
    # (A1) * X1^2 + X1^2 * (6*I) = C1  => (A1 + 6*I) * X1^2 = C1
    A1 = np.array([[5., 0.], [0., -5.]])
    C1 = np.array([[-53./12., 0.], [0., 0.]])
    I = np.identity(2)
    
    # Solve for X1_squared
    M1 = A1 + 6 * I
    inv_M1 = np.linalg.inv(M1)
    X1_squared = inv_M1 @ C1
    
    # The first coordinate of X1 is the sqrt of the (0,0) element of X1_squared
    x1_coord_sq = X1_squared[0, 0]
    sol1_x1_coord = cmath.sqrt(x1_coord_sq)
    sol2_x1_coord = -sol1_x1_coord

    # --- Second Equation ---
    # (A2) * X2^2 + X2^2 * (6*I) = C2  => (A2 + 6*I) * X2^2 = C2
    A2 = np.array([[4., 0.], [0., -5.]])
    C2 = np.array([[-3./11., 0.], [0., 0.]])
    
    # Solve for X2_squared
    M2 = A2 + 6 * I
    inv_M2 = np.linalg.inv(M2)
    X2_squared = inv_M2 @ C2
    
    # The first coordinate of X2 is the sqrt of the (0,0) element of X2_squared
    x2_coord_sq = X2_squared[0, 0]
    sol1_x2_coord = cmath.sqrt(x2_coord_sq)
    sol2_x2_coord = -sol1_x2_coord
    
    # --- Summing the coordinates ---
    total_sum = sol1_x1_coord + sol2_x1_coord + sol1_x2_coord + sol2_x2_coord
    
    print("The first equation leads to solutions for the first coordinate of X1, which are:")
    print(f"  s1 = {sol1_x1_coord}")
    print(f"  s2 = {sol2_x1_coord}")
    print("\nThe second equation leads to solutions for the first coordinate of X2, which are:")
    print(f"  s3 = {sol1_x2_coord}")
    print(f"  s4 = {sol2_x2_coord}")
    
    print("\nThe sum of the first coordinates of all solutions is:")
    print(f"Sum = ({sol1_x1_coord}) + ({sol2_x1_coord}) + ({sol1_x2_coord}) + ({sol2_x2_coord}) = {total_sum}")

solve_and_sum()
<<<0>>>