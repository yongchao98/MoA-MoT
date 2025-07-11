import math

# The theory of inhomogeneous random graphs shows that the critical time 'c'
# for the emergence of the giant component in this model is the solution
# to the following equation, where π is the mathematical constant Pi.
#
# Equation: 4 * c^2 / π^2 = 1
#
# This script solves this equation for 'c' and prints each step with the
# corresponding numerical values.

val_four = 4.0
pi = math.pi
val_one = 1.0

# Calculate π squared for use in the equation
pi_squared = pi**2

print("The equation that determines the critical time 'c' is:")
print(f"{val_four} * c^2 / π^2 = {val_one}")
print(f"\nSubstituting the numerical value for π^2 (which is approximately {pi_squared:.8f}):")
print(f"{val_four} * c^2 / {pi_squared:.8f} = {val_one}")

print("\nTo solve for 'c', we first isolate the c^2 term:")
c_squared = pi_squared / val_four
print(f"c^2 = {pi_squared:.8f} / {val_four}")
print(f"c^2 = {c_squared:.8f}")

print("\nFinally, we solve for 'c' by taking the square root:")
c = math.sqrt(c_squared)
print(f"c = √({c_squared:.8f})")
print(f"c ≈ {c:.8f}")

print("\nThe exact value of c is π/2.")
