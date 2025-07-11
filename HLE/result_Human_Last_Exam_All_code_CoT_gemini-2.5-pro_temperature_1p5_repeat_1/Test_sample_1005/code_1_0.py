import math

# Step 1-4: The analytical derivation leads to a simplified integral.
# The original integral I can be transformed through substitutions into:
# I = integral from 0 to pi/4 of (2**(-15/16) * (sin(theta))**(1/4) + 2**(17/16) * (sin(theta))**(1/4) * sec(theta)**2) d(theta)
# This can be rewritten as:
# I = integral from 0 to pi/4 of 2**(17/16) * ((sin(theta))**(1/4)) * (1/4 + sec(theta)**2) d(theta)

# Step 5: Recognize the integrand as the derivative of a function.
# The integrand is the derivative with respect to theta of the function:
# F(theta) = 2**(17/16) * tan(theta) * (sin(theta))**(1/4)
# So, the integral is F(pi/4) - F(0).

# Step 6: Final Calculation
# Define the constants and values needed for the final equation.
c = 2**(17/16)
p = 1/4
theta_upper = math.pi / 4
theta_lower = 0

print("The integral is evaluated using the Fundamental Theorem of Calculus.")
print("The antiderivative function F(theta) is C * tan(theta) * (sin(theta))^p")
print(f"where C = 2**(17/16) = {c}")
print(f"and p = 1/4 = {p}")
print(f"The integral is F(theta_upper) - F(theta_lower), with theta_upper = pi/4 and theta_lower = 0.")

# Value at the upper limit theta = pi/4
tan_upper = math.tan(theta_upper)
sin_upper = math.sin(theta_upper)
F_upper = c * tan_upper * (sin_upper)**p

print("\n--- Evaluation at the upper limit theta = pi/4 ---")
print(f"tan(pi/4) = {tan_upper}")
print(f"sin(pi/4) = {sin_upper}")
print(f"F(pi/4) = 2**(17/16) * {tan_upper} * ({sin_upper})**(1/4) = {F_upper}")

# Value at the lower limit theta = 0
tan_lower = math.tan(theta_lower)
sin_lower = math.sin(theta_lower)
# F(0) = c * tan(0) * (sin(0))**p = 0
F_lower = 0

print("\n--- Evaluation at the lower limit theta = 0 ---")
print(f"tan(0) = {tan_lower}")
# Note: although sin(0)**(1/4) involves 0**0.25 which is 0, tan(0)=0 makes the whole term 0.
print(f"F(0) = {F_lower}")

# The value of the definite integral
result = F_upper - F_lower

print("\n--- Final Result ---")
print(f"The value of the definite integral is F(pi/4) - F(0) = {result} - {F_lower} = {result}")

# The exact analytical form of the result is 2**(15/16)
exact_base = 2
exact_exponent_num = 15
exact_exponent_den = 16
print(f"\nThe exact value is {exact_base}^({exact_exponent_num}/{exact_exponent_den}).")