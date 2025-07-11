import math
import numpy as np

# Step 1 & 2: Define the characteristic quantities of the identified group Th
# R1: Sum of character table entries
R1 = 16
# R2: Number of irreducible representations
R2 = 8
# R3: Order of the group
R3 = 24
# R4: Exponent of the group
R4 = 6

print(f"The identified characteristic quantities of the group are:")
print(f"R1 (Sum of character table entries) = {R1}")
print(f"R2 (Number of irreducible representations) = {R2}")
print(f"R3 (Order of the group) = {R3}")
print(f"R4 (Exponent of the group) = {R4}")
print("-" * 30)

# Step 3: Calculate the coefficient C
# C is the floor of the contraharmonic mean of the four R values.
# Contraharmonic mean = (R1^2 + R2^2 + R3^2 + R4^2) / (R1 + R2 + R3 + R4)
C_numerator = R1**2 + R2**2 + R3**2 + R4**2
C_denominator = R1 + R2 + R3 + R4
C = math.floor(C_numerator / C_denominator)

print(f"The numerator for C is: {R1}^2 + {R2}^2 + {R3}^2 + {R4}^2 = {C_numerator}")
print(f"The denominator for C is: {R1} + {R2} + {R3} + {R4} = {C_denominator}")
print(f"C = floor({C_numerator} / {C_denominator}) = floor({C_numerator/C_denominator}) = {C}")
print("-" * 30)

# Step 4: Define polynomials Q(x) and S(x) based on C.
# Q(x) = C * (-x^2 + x^4 - x^6 + x^8)
# S(x) = C * (x - x^3 + x^5 - x^7 + x^9)
# Coefficients will be used in matrix construction.
q_coeffs = {8: C, 6: -C, 4: C, 2: -C} # deg Q = 8
s_coeffs = {9: C, 7: -C, 5: C, 3: -C, 1: C} # deg S = 9

# Step 5: Calculate the traces of M1 and M2.
# Calculate Tr(M2) for M2 = Sm(Q(x), x^10 + S(x), x)
# P1 = Q(x), n = 8. Leading coeff a_n = q_8 = C
# P2 = x^10 + S(x), m = 10. Leading coeff b_m = 1
n = 8
m = 10
a_n = q_coeffs.get(n, 0)
b_m = 1 # Coefficient of x^10 in x^10 + S(x)
trace_M2 = m * a_n + n * b_m

print("Calculating Trace(M2):")
print(f"M2 is the Sylvester Matrix of Q(x) (deg n={n}) and x^10+S(x) (deg m={m}).")
print(f"Leading coefficient of Q(x) is a_n = {a_n}")
print(f"Leading coefficient of x^10+S(x) is b_m = {b_m}")
print(f"Tr(M2) = m * a_n + n * b_m = {m} * {a_n} + {n} * {b_m} = {trace_M2}")
print("-" * 30)

# Calculate Tr(M1) for M1 = Bm(Q(x), S(x), x)
# From the derivation, Tr(Bm(Q,S)) = C^2 * (-4)
# This comes from Tr(Bm(Q/C, S/C)) = Tr(Bm(x^10-x(S/C), S/C))
# Which simplifies to -Tr(Bm(x(S/C), S/C)) = -sum_of_squares_of_coeffs(S/C)
s0_coeffs_vec = np.zeros(9)
for i in range(9):
    s0_coeffs_vec[i] = s_coeffs.get(i, 0) / C

sum_sq_s0 = np.sum(s0_coeffs_vec**2)
trace_M1 = C**2 * (-sum_sq_s0)

print("Calculating Trace(M1):")
print(f"Tr(M1) is based on the identity Q(x) + x*S(x) = {C}*x^10.")
print(f"Tr(Bm(Q/C, S/C)) equals the negative sum of squares of coefficients of S(x)/C.")
print(f"Coefficients of S(x)/C up to degree 8 are: {s0_coeffs_vec.tolist()}")
print(f"Sum of squares = {sum_sq_s0}")
print(f"Tr(M1) = C^2 * (-sum_of_squares) = {C}^2 * (-{sum_sq_s0}) = {trace_M1}")
print("-" * 30)

# Step 6: Calculate the final trace T
# T = Tr(M1 x I2 + M2) = 2 * Tr(M1) + Tr(M2)
trace_I2 = 2
T = trace_I2 * trace_M1 + trace_M2

print("Calculating the final trace T:")
print(f"T = Tr(M1 x I2) + Tr(M2)")
print(f"T = Tr(M1) * Tr(I2) + Tr(M2)")
print(f"T = {trace_M1} * {trace_I2} + {trace_M2}")
print(f"T = {trace_M1 * trace_I2} + {trace_M2}")
print(f"T = {T}")

print("\nFinal Answer:")
print(f"<<<{T}>>>")