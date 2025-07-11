# This script determines the value of a computational factor from foundational
# scientific literature on modeling phase change.

# The problem concerns a computational factor, 'C', used in the Carman-Kozeny
# source term for simulating melting. We need to find its value in a specific
# "prior, simulation-only" publication before it was modified for a later benchmark study.

# Based on a review of the literature:
# - The prior work is Voller and Prakash (1987), "A fixed grid numerical modelling methodology...".
# - The later benchmark work is Brent, Voller, and Reid (1988), "Enthalpy-porosity technique...".

# In the 1987 prior work, the authors state for their test problems, a value of 1.6 x 10^9 was used.
# In the 1988 benchmark work, this was modified to 1.6 x 10^6 to match experiments.

# The question asks for the original value from the prior (1987) paper.
# The final equation for this value is C = 1.6 * 10^9.

coefficient = 1.6
base = 10
exponent = 9

print("The question asks for the original computational factor 'C' from the prior, simulation-only publication.")
print("This value is found in the 1987 paper by Voller and Prakash.")
print("\nThe value is represented by the equation:")
# Using ** for exponentiation in Python
print(f"C = {coefficient} * {base}**{exponent}")

print("\nEach number in the final equation is:")
print(f"Coefficient: {coefficient}")
print(f"Base: {base}")
print(f"Exponent: {exponent}")

# Calculate the final value
final_value = coefficient * (base ** exponent)
print(f"\nThus, the original value used for the computational factor C was: {final_value:.1e}")
