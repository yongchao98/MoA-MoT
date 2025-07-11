# Plan:
# 1. Define the exponents of the monomials in the polynomial phase.
# 2. Choose the scaling weights based on the Newton Polyhedron.
#    The vertices are (3,0) and (0,3), suggesting weights (1,1).
# 3. Calculate the weighted degrees (which are just the total degrees here) for each monomial.
# 4. Sum these degrees.
# 5. Apply the formula for the critical exponent p.

# Step 1: Define the exponents (k_x, k_y) for each monomial
# a_1*x, a_2*y, a_3*x^2, a_4*xy, a_5*y^2, a_6*x^3, a_7*x^2*y, a_8*x*y^2, a_9*y^3
exponents = [
    (1, 0), (0, 1),  # degree 1
    (2, 0), (1, 1), (0, 2),  # degree 2
    (3, 0), (2, 1), (1, 2), (0, 3)   # degree 3
]

# Step 2: Define the scaling weights
delta_x = 1
delta_y = 1
print(f"Using scaling weights delta_x = {delta_x}, delta_y = {delta_y}")

# Step 3: Calculate the weighted degrees
nu = [k_x * delta_x + k_y * delta_y for k_x, k_y in exponents]
print("The degrees of the monomials are:")
for i, d in enumerate(nu):
    print(f"  nu_{i+1} = {d}")

# Step 4: Sum the degrees
sum_nu = sum(nu)
print(f"The sum of the degrees is: sum(nu) = {sum_nu}")

# Step 5: Apply the formula for p
sum_delta = delta_x + delta_y
print(f"The sum of the weights is: delta_x + delta_y = {sum_delta}")

p = sum_nu / sum_delta
print(f"The critical exponent p is given by the formula: p = sum(nu) / (delta_x + delta_y)")
print(f"So, p = {sum_nu} / {sum_delta} = {p}")

print("\nThe largest p such that the function I is not in L^p(R^9) is 10.")