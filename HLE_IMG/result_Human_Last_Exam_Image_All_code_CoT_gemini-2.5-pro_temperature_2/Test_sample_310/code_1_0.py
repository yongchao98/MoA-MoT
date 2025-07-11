#
# This script calculates the number-average degree of polymerization (Nn)
# for the missing simulation plot based on the provided problem description.
#

# --- Step 1: Define the parameters for the missing simulation ---

# From analyzing the plots, the initial monodisperse polymer length (N) is 20.
N = 20

# The analysis of the nine provided plots reveals that the missing plot
# corresponds to a LINEAR polymer with a degree of destruction 'a'
# associated with m=5.
m = 5
a = m / 25.0

# --- Step 2: Define the formula for Nn of a linear polymer ---

# The number-average degree of polymerization for a linear polymer undergoing
# random chain scission is given by the formula:
# Nn = N / (1 + (N - 1) * a)

# --- Step 3: Perform the calculation ---

# We plug the values of N and a into the formula.
numerator = N
denominator = 1 + (N - 1) * a
Nn_value = numerator / denominator

# --- Step 4: Print the result with the full equation ---
print("Calculation for the Number-Average Degree of Polymerization (Nn) of the Missing Plot")
print("------------------------------------------------------------------------------------")
print(f"Missing simulation parameters: Polymer Type = Linear, N = {N}, a = {a}")
print("\nFormula: Nn = N / (1 + (N - 1) * a)")
print("\nStep-by-step calculation:")
print(f"Nn = {N} / (1 + ({N} - 1) * {a})")
print(f"Nn = {N} / (1 + {N-1} * {a})")
print(f"Nn = {N} / (1 + {(N-1)*a})")
print(f"Nn = {N} / {denominator}")
print(f"\nThe final calculated value is:")
print(f"Nn = {Nn_value}")
