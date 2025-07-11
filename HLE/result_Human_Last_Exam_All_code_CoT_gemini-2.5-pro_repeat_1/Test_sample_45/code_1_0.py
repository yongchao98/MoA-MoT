import sympy

# Define symbols
f, t, n, g, omega, P, I = sympy.symbols('f t n g omega P I')

print("This script derives the exponent 'n' in the relationship f is proportional to t^n for a liquid-mirror telescope spun up by a constant power source.")
print("-" * 80)

# Step 1: Relate focal length (f) to angular velocity (omega)
print("Step 1: The focal length (f) of a rotating liquid mirror is given by f = g / (2 * omega^2).")
print("This means f is proportional to omega^(-2).")
print("-" * 80)

# Step 2: Relate angular velocity (omega) to time (t)
print("Step 2: A constant power source (P) increases the rotational kinetic energy (E_k) such that E_k = P * t.")
print("The kinetic energy is also given by E_k = (1/2) * I * omega^2.")
print("Equating these gives (1/2) * I * omega^2 = P * t.")
print("Solving for omega^2, we find omega^2 = (2 * P * t) / I.")
print("This means omega^2 is proportional to t.")
print("-" * 80)

# Step 3: Combine the relationships to find f(t) and solve for n
print("Step 3: We substitute the time-dependence of omega^2 into the equation for f.")
print("f is proportional to omega^(-2), which is proportional to (t)^(-1).")
print("So, f is proportional to t^(-1).")
print("-" * 80)

# Final conclusion
print("We are looking for n where f is proportional to t^n.")
final_n = -1
print(f"Comparing f proportional to t^({final_n}) with f proportional to t^n, we find:")
print(f"n = {final_n}")