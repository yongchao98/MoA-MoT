import numpy as np
from sympy import S, dedekind_sum

# Step 1: Define the matrices for the Dehn twists
Ma = np.array([[1, 1], [0, 1]])
Mb = np.array([[1, 0], [1, 1]])

# Step 2: Compute the matrix for the element psi = (Da Db)^3
C = np.dot(Ma, Mb)
M_psi = np.linalg.matrix_power(C, 3)

p, q = M_psi[0, 0], M_psi[0, 1]
r, s = M_psi[1, 0], M_psi[1, 1]

print(f"The matrix for psi = (Da * Db)^3 is:\n{M_psi}")
print("-" * 30)

# Step 3: Extract components and calculate the trace
trace_M = np.trace(M_psi)
sign_r = np.sign(r)

# Calculate the Dedekind sum s(s, |r|)
# Note: Dedekind sums are defined for coprime integers.
# The entries of powers of this matrix have this property.
dedekind_val = dedekind_sum(s, r)

print("Calculating the FDTC for psi = (Da * Db)^3 using the formula:")
print("tau(psi) = (1/12) * (tr(M) - 3*sign(r) - 12*s(s, |r|))")
print(f"  tr(M) = {trace_M}")
print(f"  r = {r}, sign(r) = {int(sign_r)}")
print(f"  s = {s}")
print(f"  Dedekind sum s({s}, {r}) = {dedekind_val}")
print("-" * 30)

# Step 4: Calculate the FDTC for psi
tau_psi = (S(1)/12) * (trace_M - 3 * sign_r - 12 * dedekind_val)

print(f"The Fractional Dehn Twist Coefficient for (Da * Db)^3 is:")
print(f"tau(psi) = (1/12) * ({trace_M} - 3*({int(sign_r)}) - 12*({dedekind_val}))")
print(f"tau(psi) = {tau_psi}")
print("-" * 30)


# Step 5: The FDTC for phi = psi^3 is 3 * tau_psi
total_tau = 3 * tau_psi

print(f"The final FDTC for (Da * Db)^9 is 3 * tau(psi):")
print(f"tau((Da * Db)^9) = 3 * {tau_psi} = {total_tau}")
