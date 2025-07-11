import numpy as np

# This function solves for the first coordinates and returns them
def solve_first_coordinates(A, C):
    """
    Solves the matrix equation (A + 6I)X^2 = C for the first coordinates of X.
    """
    # The matrix B is 6*I, where I is the 2x2 identity matrix
    I = np.identity(2)
    
    # Simplify the equation to M * X_sq = C, where M = A + 6I
    M = A + 6 * I
    
    # Solve for X_sq = inv(M) * C
    inv_M = np.linalg.inv(M)
    X_sq = inv_M @ C
    
    # The resulting X_sq is of the form [[p, 0], [0, 0]].
    # For X^2 = X_sq, where X = [[a, b], [c, d]], we must have b=c=d=0.
    # This leads to a^2 = p, where p is the top-left element of X_sq.
    p = X_sq[0, 0]
    
    # The solutions for the first coordinate 'a' are sqrt(p) and -sqrt(p).
    # Since p is negative, we use complex numbers.
    sol1 = np.sqrt(complex(p))
    sol2 = -sol1
    
    return sol1, sol2

# --- Equation 1 ---
A1 = np.array([[5, 0], [0, -5]])
C1 = np.array([[-53/12, 0], [0, 0]])
sol_x1_1, sol_x1_2 = solve_first_coordinates(A1, C1)

# --- Equation 2 ---
A2 = np.array([[4, 0], [0, -5]])
C2 = np.array([[-3/11, 0], [0, 0]])
sol_x2_1, sol_x2_2 = solve_first_coordinates(A2, C2)

# The total sum is the sum of all possible first coordinates
total_sum = sol_x1_1 + sol_x1_2 + sol_x2_1 + sol_x2_2

# Print the final equation with all the numbers
print("The first coordinates for X1 are: {} and {}".format(sol_x1_1, sol_x1_2))
print("The first coordinates for X2 are: {} and {}".format(sol_x2_1, sol_x2_2))
print("\nThe final equation for the sum is:")
print("({:.8f}) + ({:.8f}) + ({:.8f}) + ({:.8f}) = {:.8f}".format(sol_x1_1, sol_x1_2, sol_x2_1, sol_x2_2, total_sum))
