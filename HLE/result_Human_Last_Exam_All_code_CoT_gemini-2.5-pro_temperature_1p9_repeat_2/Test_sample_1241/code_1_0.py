import numpy as np

# Define the given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

print("Solving for the steady-state probabilities p0, p1, p2, p3.")
print("The system of equations is:")
print(f"1) 0 = -{lambda_01}*p0 + {lambda_10}*p1")
print(f"2) 0 = {lambda_12}*p1 - ({lambda_21} + {lambda_23})*p2")
print(f"3) 0 = {lambda_23}*p2 - {lambda_31}*p3")
print("4) 1 = p0 + p1 + p2 + p3")
print("-" * 30)

# Step 1: Express p0, p2, p3 in terms of p1
# From equation 1: p0 = (lambda_10 / lambda_01) * p1
k0 = lambda_10 / lambda_01
print(f"From eq 1, p0 = ({lambda_10} / {lambda_01}) * p1 = {k0:.4f} * p1")

# From equation 2: p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
k2_num = lambda_12
k2_den = lambda_21 + lambda_23
k2 = k2_num / k2_den
print(f"From eq 2, p2 = ({lambda_12} / ({lambda_21} + {lambda_23})) * p1 = {k2:.4f} * p1")

# From equation 3: p3 = (lambda_23 / lambda_31) * p2
k3 = lambda_23 / lambda_31
# Since p2 = k2 * p1, then p3 = k3 * k2 * p1
print(f"From eq 3, p3 = ({lambda_23} / {lambda_31}) * p2 = {k3:.4f} * p2 = {k3*k2:.4f} * p1")
print("-" * 30)

# Step 2: Substitute into the normalization equation to find p1
# p0 + p1 + p2 + p3 = 1
# (k0 * p1) + p1 + (k2 * p1) + (k3 * k2 * p1) = 1
# p1 * (k0 + 1 + k2 + k3 * k2) = 1
p1_coeff = k0 + 1 + k2 + (k3 * k2)
p1 = 1 / p1_coeff

# Step 3: Calculate p0
p0 = k0 * p1

# Step 4: Calculate the required sum p0 + p1
result = p0 + p1

print("Solving for p1 using the normalization equation p0 + p1 + p2 + p3 = 1:")
print(f"p1 * ({k0:.4f} + 1 + {k2:.4f} + {k3*k2:.4f}) = 1")
print(f"p1 * {p1_coeff:.4f} = 1")
print(f"p1 = {p1:.10f}")
print("-" * 30)

print(f"Now we calculate p0:")
print(f"p0 = {k0:.4f} * p1 = {k0:.4f} * {p1:.10f} = {p0:.10f}")
print("-" * 30)

print("Finally, we calculate the required sum P0(inf) + P1(inf), which is p0 + p1.")
print("The final equation is:")
print(f"p0 + p1 = {p0:.10f} + {p1:.10f}")
print(f"p0 + p1 = {result}")
