import scipy.constants

# The user wants to find the ratio of the uncertainty of the momentum (delta_p)
# of an electron to its momentum (p) in the first Bohr orbit.

# 1. Define the given and known constants in picometers (pm) to keep units consistent.
# The uncertainty of the electron's position is given. 1 pm = 1e-12 m.
delta_x_pm = 10.0

# The radius of the first Bohr orbit (n=1) is the Bohr radius, a0.
# We get this value from scipy.constants and convert it to pm.
# scipy.constants.physical_constants['Bohr radius'][0] gives the value in meters.
a0_m = scipy.constants.physical_constants['Bohr radius'][0]
a0_pm = a0_m * 1e12

# 2. Formulate the ratio.
# From the Heisenberg Uncertainty Principle, the minimum uncertainty in momentum is:
# delta_p = hbar / (2 * delta_x)
# From the Bohr model, the momentum of an electron in the first orbit (n=1) is:
# p = hbar / a0 (where a0 is the radius of the first orbit)
#
# The ratio is delta_p / p = (hbar / (2 * delta_x)) / (hbar / a0)
# The hbar term cancels out, simplifying the ratio to: a0 / (2 * delta_x)

# 3. Calculate the ratio using the values in picometers.
ratio = a0_pm / (2 * delta_x_pm)

# 4. Print the final equation and the result.
print("This script calculates the ratio of the uncertainty of an electron's momentum to its momentum in the first Bohr orbit.")
print("The ratio is determined by the formula: a0 / (2 * Δx)")
print("where a0 is the Bohr radius and Δx is the uncertainty in position.\n")

print("Equation setup:")
print(f"Ratio = (Bohr Radius) / (2 * Uncertainty in Position)")
print(f"Ratio = {a0_pm:.4f} pm / (2 * {delta_x_pm} pm)")
print(f"Ratio = {a0_pm:.4f} / {2 * delta_x_pm}")
print(f"The calculated ratio is: {ratio:.4f}")