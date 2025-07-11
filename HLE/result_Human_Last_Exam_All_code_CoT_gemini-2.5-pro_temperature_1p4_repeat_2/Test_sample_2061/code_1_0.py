import math

# Given constants
T = math.log(10)
B = 0.5 * (10**20) / (99**2)

# Calculate intermediate constants k2 and k3
q = math.exp(-T) # q = 0.1
k2 = 1 - q**2    # 1 - 0.01 = 0.99
k3 = 1 - q**3    # 1 - 0.001 = 0.999

# Numerator of the expression for alpha^2
numerator = 2 * B * (k2**2) * k3

# Denominator of the expression for alpha^2
denominator = 3 * k2 - k3

# Calculate alpha^2 and then alpha
alpha_squared = numerator / denominator
alpha = math.sqrt(alpha_squared)

print(f"T = {T}")
print(f"B = {B}")
print(f"e^-T = {q}")
print(f"k2 = 1 - e^(-2T) = {k2}")
print(f"k3 = 1 - e^(-3T) = {k3}")
print(f"Numerator = 2 * B * k2^2 * k3 = {numerator}")
print(f"Denominator = 3*k2 - k3 = {denominator}")
print(f"alpha^2 = Numerator / Denominator = {alpha_squared}")
print(f"alpha = sqrt(alpha^2) = {alpha}")

# Final equation demonstration
# The value of alpha squared is exactly 37/73 * 10^16
alpha_sq_exact_num = 37 * 10**16
alpha_sq_exact_den = 73
alpha_exact = 10**8 * math.sqrt(37/73)
print(f"The exact value of alpha is 10^8 * sqrt(37/73)")
print(f"alpha = {alpha_exact}")

# To satisfy the output format, we must present the final equation
# The equation for alpha is alpha = sqrt( (2*B*k2^2*k3)/(3*k2-k3) )
print("The final equation for calculating alpha is: ")
print(f"alpha = sqrt( (2 * {B} * {k2}^2 * {k3}) / (3 * {k2} - {k3}) )")
print(f"alpha = {alpha}")