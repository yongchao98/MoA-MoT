# The final analytical expression is derived from the critical-state model
# for a thin rectangular superconductor in a perpendicular magnetic field.

# The derived formula is:
# M(H) = - (Jc * a / 2) * (1 - exp(- (pi * H) / (b * Jc)))

print("The analytical expression for the initial magnetization curve M(H) is:")
# We print the formula as a string.
print("M(H) = - (Jc * a / 2) * (1 - exp(- (pi * H) / (b * Jc)))")

print("\n--- Equation Breakdown ---")

# As requested, here is a breakdown of each part of the final equation,
# including the numerical constants.
print("1. The overall sign is negative, indicated by the constant '-1'. This reflects the diamagnetic nature of the superconductor's response.")
print("\n2. The saturation magnetization term is 'Jc * a / 2'.")
print("   - This term represents the maximum absolute value of magnetization.")
print("   - It involves the numerical constant '2' in the denominator.")

print("\n3. The approach-to-saturation term is '(1 - exp(...))'.")
print("   - This term starts at 0 (for H=0) and goes to 1 as H increases.")
print("   - It involves the numerical constant '1'.")

print("\n4. The argument of the exponential function 'exp' is '- (pi * H) / (b * Jc)'.")
print("   - This term involves the mathematical constant 'pi' (approx. 3.14159).")

print("\n--- Variable Definitions ---")
print("  M: Magnetization")
print("  H: Applied magnetic field")
print("  Jc: Critical current density")
print("  a: Half-width of the superconducting strip (long side)")
print("  b: Half-thickness of the superconducting strip (short side)")