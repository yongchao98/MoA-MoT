import numpy as np
import cmath

def solve_and_sum():
    """
    Solves the given system of matrix equations and finds the sum of the
    first coordinates of all solutions.
    """
    # --- Equation 1 ---
    # (5 0; 0 -5) X1^2 + X1^2 (6 0; 0 6) = (-53/12 0; 0 0)
    # This is an equation of the form A1*Y + Y*B = C1, where Y = X1^2.
    # Since B = 6*I, this simplifies to (A1 + 6*I)*Y = C1.
    A1 = np.array([[5., 0.], [0., -5.]])
    C1 = np.array([[-53./12., 0.], [0., 0.]])
    I = np.identity(2)

    # Solve for Y = X1^2
    M1 = A1 + 6 * I
    M1_inv = np.linalg.inv(M1)
    X1_sq = M1_inv @ C1

    # The first coordinate of X1 is a, where a^2 = X1_sq[0, 0].
    p1 = X1_sq[0, 0]
    sol1_a_pos = cmath.sqrt(p1)
    sol1_a_neg = -sol1_a_pos

    # --- Equation 2 ---
    # (4 0; 0 -5) X2^2 + X2^2 (6 0; 0 6) = (-3/11 0; 0 0)
    # Similar simplification applies.
    A2 = np.array([[4., 0.], [0., -5.]])
    C2 = np.array([[-3./11., 0.], [0., 0.]])

    # Solve for Y = X2^2
    M2 = A2 + 6 * I
    M2_inv = np.linalg.inv(M2)
    X2_sq = M2_inv @ C2
    
    # The first coordinate of X2 is a, where a^2 = X2_sq[0, 0].
    p2 = X2_sq[0, 0]
    sol2_a_pos = cmath.sqrt(p2)
    sol2_a_neg = -sol2_a_pos

    # The solutions for the first coordinate of X1 are the two square roots of p1.
    print(f"The two solutions for the first coordinate of X_1 are: {sol1_a_pos} and {sol1_a_neg}")
    # The solutions for the first coordinate of X2 are the two square roots of p2.
    print(f"The two solutions for the first coordinate of X_2 are: {sol2_a_pos} and {sol2_a_neg}")
    
    # The total sum is the sum of all four solution coordinates.
    total_sum = sol1_a_pos + sol1_a_neg + sol2_a_pos + sol2_a_neg
    
    print("\nThe sum of the first coordinate of all solutions is calculated as follows:")
    # Using format to make it clear what is being added.
    sum_eq = f"({sol1_a_pos}) + ({sol1_a_neg}) + ({sol2_a_pos}) + ({sol2_a_neg})"
    print(f"{sum_eq} = {total_sum}")

solve_and_sum()