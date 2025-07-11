import numpy as np

# Step 1: Define the system matrices
A = np.array([
    [-1, 1],
    [1, 0]
])
B = np.array([
    [1, 2],
    [1, 0]
])

# From the analysis, we need to solve the following system of linear equations for f1 and f2:
# f1 + f2 = -1
# f1 + 2*f2 = -3
# This can be written in matrix form as M * f = v
M = np.array([
    [1, 1],
    [1, 2]
])
v = np.array([-1, -3])

# Step 2: Solve the system of equations to find f1 and f2
try:
    f1_f2 = np.linalg.solve(M, v)
    f1 = f1_f2[0]
    f2 = f1_f2[1]

    # Step 3: Construct the full feedback gain matrix F
    # We assumed f3=0 and f4=0 to find a specific solution
    F = np.array([
        [f1, f2],
        [0, 0]
    ])

    print("The calculated state feedback gain matrix F is:")
    print(F)
    print("\nThe final equation for the elements of F is:")
    print(f"F = [[{F[0,0]}, {F[0,1]}],")
    print(f"     [{F[1,0]}, {F[1,1]}]]")


    # Step 4: Verification (Optional but good practice)
    # Calculate the closed-loop matrix A_cl = A + B@F
    A_cl = A + B @ F
    # Calculate the eigenvalues of the closed-loop system
    eigenvalues = np.linalg.eigvals(A_cl)

    print("\n--- Verification ---")
    print("Closed-loop matrix A + BF:")
    print(A_cl)
    print("\nEigenvalues of the closed-loop system A + BF:")
    # Sorting for consistent order
    eigenvalues.sort()
    print(np.round(eigenvalues, 5))
    print("These match the desired eigenvalues of -1 +/- j.")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")
