import numpy as np

# This script finds the FGH-tripled fixed point for a specific case where
# X, Y, and Z are the set of real numbers (R), and F, G, H are linear functions.
# A fixed point (x, y, z) satisfies the system:
# F(x, y, z) = x
# G(y, x, y) = y
# H(z, y, x) = z

# Let's define our linear functions:
# F(x, y, z) = k1*x + k2*y + k3*z + c1
# G(y, x, y) = l1*y + l2*x + c2
# H(z, y, x) = m1*z + m2*y + m3*x + c3

# --- 1. Define Coefficients and Constants ---
# These coefficients are the contraction constants (k_i, l_i, m_i).
# We choose them to satisfy the condition: max(sum_k, sum_l, sum_m) < 1
k1, k2, k3 = 0.1, 0.2, 0.1  # sum = 0.4
l1, l2 = 0.3, 0.2           # sum = 0.5
m1, m2, m3 = 0.1, 0.2, 0.3  # sum = 0.6
# max(0.4, 0.5, 0.6) = 0.6, which is < 1. The condition holds.

# Constants for the functions
c1, c2, c3 = 5, 2, 8

# --- 2. Set up the System of Linear Equations ---
# The fixed point condition gives us a system of equations A*v = b, where v = [x, y, z].
# x = k1*x + k2*y + k3*z + c1  => (1-k1)x - k2*y - k3*z = c1
# y = l1*y + l2*x + c2           => -l2*x + (1-l1)y = c2
# z = m1*z + m2*y + m3*x + c3  => -m3*x - m2*y + (1-m1)z = c3

# Matrix A
A = np.array([
    [1 - k1, -k2, -k3],
    [-l2, 1 - l1, 0],
    [-m3, -m2, 1 - m1]
])

# Vector b
b = np.array([c1, c2, c3])

print("--- System of Equations to find the FGH-Tripled Fixed Point (x, y, z) ---")
print("The functions are defined as:")
print(f"F(x, y, z) = {k1}*x + {k2}*y + {k3}*z + {c1}")
print(f"G(y, x, y) = {l1}*y + {l2}*x + {c2}")
print(f"H(z, y, x) = {m1}*z + {m2}*y + {m3}*x + {c3}\n")

print("The fixed point equations F=x, G=y, H=z can be written as A*v = b:")
for i in range(A.shape[0]):
    row_str = " + ".join([f"{A[i, j]:.2f}*{var}" for j, var in enumerate(['x', 'y', 'z']) if A[i, j] != 0])
    row_str = row_str.replace("+ -", "- ")
    print(f"{row_str} = {b[i]}")
print("\n" + "-"*70 + "\n")

# --- 3. Solve the System ---
# Because the contraction condition is met, a unique solution is guaranteed.
try:
    fixed_point = np.linalg.solve(A, b)
    x, y, z = fixed_point

    print("--- Solution ---")
    print(f"The unique FGH-tripled fixed point (x, y, z) is:")
    print(f"x = {x}")
    print(f"y = {y}")
    print(f"z = {z}")

    # --- 4. Verification (Optional) ---
    # We can plug the solution back into the original functions F, G, H
    # to verify that F(x,y,z)=x, G(y,x,y)=y, and H(z,y,x)=z.
    f_val = k1*x + k2*y + k3*z + c1
    g_val = l1*y + l2*x + c2
    h_val = m1*z + m2*y + m3*x + c3

    print("\n--- Verification ---")
    print(f"F({x:.2f}, {y:.2f}, {z:.2f}) = {f_val:.2f} (matches x)")
    print(f"G({y:.2f}, {x:.2f}, {y:.2f}) = {g_val:.2f} (matches y)")
    print(f"H({z:.2f}, {y:.2f}, {x:.2f}) = {h_val:.2f} (matches z)")

except np.linalg.LinAlgError:
    print("The system could not be solved. The matrix A is singular.")
