import math

# Step 1 & 2: Define constants and roots of the characteristic equation
sqrt5 = math.sqrt(5)
r1 = (1 + sqrt5) / 2
r2 = (1 - sqrt5) / 2

# Step 3: Apply boundary conditions to find the constant C1
# We use t = ln(5) for the second boundary condition
# 5^r1 = exp(r1*ln(5)) and 5^r2 = exp(r2*ln(5))
t_bc = math.log(5)
term1_bc = math.exp(r1 * t_bc)  # This is 5^r1
term2_bc = math.exp(r2 * t_bc)  # This is 5^r2

# Calculate C1 using the second boundary condition
# C1 = (1/sqrt5) * (2 - (5^r1 + 5^r2)) / (5^r1 - 5^r2)
C1_numerator = 2 - (term1_bc + term2_bc)
C1_denominator = term1_bc - term2_bc
C1 = (1 / sqrt5) * (C1_numerator / C1_denominator)

# Step 4: Define the solution function phi_0(t)
# phi_0(t) = C1 * (exp(r1*t) - exp(r2*t))
def phi_0(t):
    term1_t = math.exp(r1 * t)
    term2_t = math.exp(r2 * t)
    return C1 * (term1_t - term2_t)

# Step 5: Calculate the final value
t_final = math.log(10**10)

phi_val = phi_0(t_final)
const_val = 2 / sqrt5

# The final equation is -phi_0(ln(10^10)) + 2/sqrt(5)
# We print the numbers involved
print(f"The value of -phi_0(ln(10^10)) is: {-phi_val}")
print(f"The value of 2/sqrt(5) is: {const_val}")

result = -phi_val + const_val

# Final equation with calculated values:
print(f"Final equation: ({-phi_val}) + ({const_val})")
print(f"Result: {result}")