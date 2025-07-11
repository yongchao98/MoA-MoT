import numpy as np

def solve_and_sum_first_coordinate(A, B, C, eq_num):
    """
    Solves the matrix equation A*X^2 + X^2*B = C and returns the sum
    of the possible values for the first coordinate of X.
    """
    print(f"--- Solving Equation {eq_num} ---")
    
    # B is a scalar matrix, B = k*I. Here k=6.
    # The equation simplifies to (A + k*I) * X^2 = C
    k = B[0, 0]
    M = A + k * np.identity(A.shape[0])
    
    # Solve for Y = X^2. Y = inv(M) * C
    inv_M = np.linalg.inv(M)
    Y = np.dot(inv_M, C)
    
    # Let X = [[a, b], [c, d]]. Then X^2 = [[a^2+bc, ...], ...].
    # The structure of the solved Y matrix is diagonal.
    # X^2 = Y implies that either b=0 and c=0, or a+d=0.
    # The case a+d=0 leads to a contradiction.
    # Thus, b=c=0, which means X is diagonal.
    # Therefore, a^2 = Y[0, 0].
    a_squared = Y[0, 0]
    
    print(f"The equation for the first coordinate 'a' is: a^2 = {a_squared:.4f}")
    
    # The solutions for 'a' are sqrt(a_squared) and -sqrt(a_squared).
    # We use np.sqrt which can handle complex numbers.
    sol1 = np.sqrt(a_squared, dtype=np.complex128)
    sol2 = -sol1
    
    print(f"The two solutions for the first coordinate are: {sol1} and {sol2}")
    
    # The sum of these two solutions is 0.
    coord_sum = sol1 + sol2
    print(f"The sum of these solutions is: ({sol1}) + ({sol2}) = {np.real_if_close(coord_sum)}")
    
    return coord_sum

# --- Equation 1 ---
A1 = np.array([[5, 0], [0, -5]])
B1 = np.array([[6, 0], [0, 6]])
C1 = np.array([[-53/12, 0], [0, 0]])
sum1 = solve_and_sum_first_coordinate(A1, B1, C1, 1)

# --- Equation 2 ---
A2 = np.array([[4, 0], [0, -5]])
B2 = np.array([[6, 0], [0, 6]])
C2 = np.array([[-3/11, 0], [0, 0]])
sum2 = solve_and_sum_first_coordinate(A2, B2, C2, 2)

# --- Final Sum ---
total_sum = sum1 + sum2
print("\n--- Final Calculation ---")
print("The sum of the first coordinate of all solutions is the sum of the individual sums.")
# The final equation as requested: sum1 + sum2 = total_sum
print(f"Final Sum = {np.real_if_close(sum1)} + {np.real_if_close(sum2)} = {np.real_if_close(total_sum)}")

<<<0>>>