import math

# Parameters for the complete intersection X
# X is an intersection of k hypersurfaces of degrees d_i in CP^N
# d_1 = 2, d_2 = 2
degrees = [2, 2]
# N = 102
N = 102

k = len(degrees)
# Dimension of X, n = N - k
n = N - k

# --- Step 1: Calculate the degree of X ---
deg_X = math.prod(degrees)

# --- Step 2: Calculate the coefficient C = [z^n] (1+z)^(N+1) / prod(1+d_i*z) ---
# This is given by the sum:
# C = sum_{j=0 to n} C(N+1, n-j) * (-1)^j * (j+1) * 2^j
# for the specific case where d_1=d_2=2
# where n=100, N=102
# The sum arises from the product of the series for (1+z)^103 and (1+2z)^-2
C = 0
for j in range(n + 1):
    # Binomial coefficient C(n,k) is math.comb(n,k)
    term = math.comb(N + 1, n - j) * ((-1)**j) * (j + 1) * (2**j)
    C += term

# --- Step 3: Calculate the Euler characteristic chi(X) = deg(X) * C ---
chi_X = deg_X * C

# --- Step 4: Calculate the dimension of the middle cohomology group ---
# b_n(X) = chi(X) - n, for n even. Here n = 100.
b_n = chi_X - n

print(f"The complex dimension of the variety X is n = {n}.")
print(f"The degree of X is {deg_X}.")
print(f"The relevant coefficient C from the generating function is {C}.")
print(f"The Euler characteristic of X is chi(X) = {deg_X} * {C} = {chi_X}.")
print(f"The dimension of the middle cohomology group H^{n}(X, Q) is b_n(X) = chi(X) - n.")
print(f"Dimension = {chi_X} - {n} = {b_n}.")
