import math

# Step 1: Define the coupling constants for neutrinos
c_V = 0.5
c_A = 0.5

# Calculate the common coupling term (c_V^2 + c_A^2)
coupling_term = c_V**2 + c_A**2

# Step 2: Calculate X_1
# From the relationship |M|^2 = 2 * sqrt(2) * G_F * m_Z^4 * (c_V^2 + c_A^2)
# and |M|^2 = X_1 * G_F * m_Z^4, we get:
# X_1 = 2 * sqrt(2) * (c_V^2 + c_A^2)
X1 = 2 * math.sqrt(2) * coupling_term

# Step 3: Calculate X_2
# From the relationship Gamma = (G_F * m_Z^3) / (12 * sqrt(2) * pi) * (c_V^2 + c_A^2)
# and Gamma = X_2 * G_F * m_Z^3, we get:
# X_2 = (c_V^2 + c_A^2) / (12 * sqrt(2) * pi)
X2 = coupling_term / (12 * math.sqrt(2) * math.pi)

# Step 4: Calculate the final result (X_1 * X_2)^(-1)
product_X1_X2 = X1 * X2
result = 1 / product_X1_X2

# Step 5: Print the components of the final equation and the result
print(f"For the Z to neutrino decay, with c_V = {c_V} and c_A = {c_A}:")
print(f"The equation for the squared amplitude is |M|^2 = X_1 * G_F * m_Z^4")
print(f"The calculated value for X_1 is: {X1:.4f}")
print("-" * 30)
print(f"The equation for the decay rate is Gamma = X_2 * G_F * m_Z^3")
print(f"The calculated value for X_2 is: {X2:.4f}")
print("-" * 30)
print("The final calculation is (X_1 * X_2)^(-1)")
print(f"Result = ({X1:.4f} * {X2:.4f})^(-1) = ({product_X1_X2:.4f})^(-1)")
print(f"The final numerical result is: {result:.4f}")

# The exact result is 24 * pi
exact_result = 24 * math.pi
print(f"The exact symbolic result is 24 * pi, which is approximately: {exact_result:.4f}")

<<<75.3982>>>