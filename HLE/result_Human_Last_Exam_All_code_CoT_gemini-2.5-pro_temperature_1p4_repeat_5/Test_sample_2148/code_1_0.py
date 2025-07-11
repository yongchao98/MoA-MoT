import math

# Step 1: Define the given coupling constants for neutrinos
cV = 0.5
cA = 0.5
c_sq_sum = cV**2 + cA**2
print(f"Given coupling constants: c_V = {cV}, c_A = {cA}")
print(f"The sum of their squares is: c_V^2 + c_A^2 = {c_sq_sum}\n")

# Step 2: Determine X1 by equating two expressions for the decay rate Gamma.
# Expression 1 (from phase space and the definition of X1): Gamma = (X1 / (16*pi)) * G_F * m_Z^3
# Expression 2 (given in the problem): Gamma = (G_F*m_Z^3 / (12*sqrt(2)*pi)) * (c_V^2 + c_A^2)
# Equating them gives: X1 / (16*pi) = (c_V^2 + c_A^2) / (12*sqrt(2)*pi)
# Solving for X1: X1 = 16 * (c_V^2 + c_A^2) / (12*sqrt(2))
X1 = (16 * c_sq_sum) / (12 * math.sqrt(2))
print(f"Calculation of X1:")
print(f"X1 = 16 * ({c_sq_sum}) / (12 * sqrt(2))")
print(f"Value of X1 = {X1}\n")


# Step 3: Determine X2 from its definition.
# Definition: Gamma = X2 * G_F * m_Z^3
# Comparing this with Expression 2 for Gamma:
# X2 = (c_V^2 + c_A^2) / (12*sqrt(2)*pi)
X2 = c_sq_sum / (12 * math.sqrt(2) * math.pi)
print(f"Calculation of X2:")
print(f"X2 = ({c_sq_sum}) / (12 * sqrt(2) * pi)")
print(f"Value of X2 = {X2}\n")


# Step 4: Calculate the product X1 * X2
X1_X2 = X1 * X2
print(f"Calculation of the product X1 * X2:")
print(f"X1 * X2 = ({X1}) * ({X2})")
print(f"Value of X1 * X2 = {X1_X2}\n")


# Step 5: Calculate the inverse (X1 * X2)^-1
# The symbolic calculation is:
# X1 * X2 = [16*0.5 / (12*sqrt(2))] * [0.5 / (12*sqrt(2)*pi)]
#         = [8 / (12*sqrt(2))] * [0.5 / (12*sqrt(2)*pi)]
#         = [2 / (3*sqrt(2))] * [1 / (24*sqrt(2)*pi)]
#         = 2 / (3*24*2*pi) = 1 / (72*pi)
# So, (X1 * X2)^-1 = 72 * pi.
inverse_X1_X2 = 1 / X1_X2
final_symbolic_value = "72 * pi"
print(f"Final calculation of (X1 * X2)^-1:")
print(f"The final equation is: (X1 * X2)^-1 = ({X1_X2})^-1 = {inverse_X1_X2}")
print(f"This corresponds to the symbolic value: {final_symbolic_value}")