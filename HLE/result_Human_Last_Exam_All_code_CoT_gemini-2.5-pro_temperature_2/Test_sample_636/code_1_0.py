import numpy as np

# Define the system matrices based on the problem statement
A = np.array([[-1, 1], [1, 0]])
B = np.array([[1, 2], [1, 0]])

# Define the calculated state feedback gain matrix F
# From our derivation, a valid F is [[1, -2], [0, 0]]
F = np.array([[1.0, -2.0], [0.0, 0.0]])

# Calculate the closed-loop matrix A_cl = A + B*F
A_cl = A + B @ F

# Calculate the eigenvalues of the closed-loop system to verify
eigenvalues = np.linalg.eigvals(A_cl)

print("Based on the derivation, a valid state feedback gain matrix F is:")
print(F)
print("\nThe final equation A + B*F = A_cl is:")

# Use array2string for clean, aligned matrix printing
print(f"{np.array2string(A, sign=' ')} \n\n+ \n\n{np.array2string(B, sign=' ')} \n\n* \n\n{np.array2string(F, sign=' ')} \n\n= \n\n{np.array2string(A_cl, sign=' ')}")


print("\nThe eigenvalues of the resulting closed-loop system (A + BF) are:")
# Format the complex numbers for clean output
# Using np.round to handle potential floating point inaccuracies
rounded_eigenvalues = np.round(eigenvalues, 5)
formatted_eigenvalues = [f"{v.real}{v.imag:+.1f}j" if v.imag != 0 else f"{v.real}" for v in rounded_eigenvalues]
print(formatted_eigenvalues)

print("\nThe eigenvalues match the desired poles at -1 + j and -1 - j.")