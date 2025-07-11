import numpy as np

# Step 1: Solve for the coefficients of Phi(z) = a*z + b
# The system of linear equations is:
# 8*a + 1*b = 5
# 2*a + 0*b = 1

# Matrix form: A * x = B
# where x = [a, b]
A = np.array([[8, 1], [2, 0]])
B = np.array([5, 1])

# Solve for a and b
try:
    coeffs = np.linalg.solve(A, B)
    a = coeffs[0]
    b = coeffs[1]
    print(f"Solving the system of equations for the coefficients of Phi(z):")
    print(f"8a + b = 5")
    print(f"2a = 1")
    print(f"The solution is a = {a}, b = {b}")
    print(f"Therefore, Phi(z) = {a}*z + {b}\n")
except np.linalg.LinAlgError:
    print("The system of equations could not be solved.")
    a, b = 0, 0 # Fallback

# Step 2: Calculate the minimum total energy
# E_total = integral from 0 to 1 of (Phi(z(x))) dx
# where z(x) = T1(sqrt(2)*x) + T2(x)
# The integrand is Phi(z) = a*z + b.
# For non-negative energy density, Phi(z) >= 0 => a*z + b >= 0.
# 0.5*z + 1 >= 0 => z >= -2.
# The minimum value of z(x) is -2.
z_min = -2

# The minimum value of the integrand is when z(x) is at its minimum.
min_integrand_value = a * z_min + b

# The minimum total energy is the integral of this minimum integrand value.
# E_total_min = integral from 0 to 1 of (min_integrand_value) dx
# This is simply min_integrand_value * (1 - 0)
E_total_min = min_integrand_value * 1

print("The minimum possible value of z(x) = T1(sqrt(2)*x) + T2(x) is z_min = -2.")
print(f"The minimum value of the integrand Phi(z_min) is {a} * ({z_min}) + {b} = {min_integrand_value}")
print(f"The minimum total energy is the integral of this constant value from 0 to 1.")
print(f"E_total_min = {min_integrand_value} * (1 - 0)")
print(f"\nThe minimum value of the total heat energy is {E_total_min}.")
