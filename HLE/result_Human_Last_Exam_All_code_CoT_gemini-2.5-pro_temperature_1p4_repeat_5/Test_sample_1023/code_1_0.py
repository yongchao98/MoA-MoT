import sys

# Step 1: Define the fundamental evolutionary rates.
# We are calculating a ratio, so we can set mu to a nominal value of 1.0.
# All other rates will be determined relative to mu.
mu = 1.0

# Step 2: Apply the assumptions from the problem to relate the rates.
# Assumption 1: The system is in equilibrium (origination rate = extinction rate).
# This implies that the speciation rate (lambda) equals the extinction rate (mu).
lambda_rate = mu

# Assumption 2: A plausible interpretation of "all processes ... occur at the same rates"
# is that the evolutionary turnover rate (lambda + mu) equals the total morphospecies
# extinction rate (mu + phi + 0.5*lambda).
# (lambda + mu) = (mu + phi + 0.5*lambda)  which simplifies to ->  phi = 0.5 * lambda.
phi_rate = 0.5 * lambda_rate

# Step 3: Calculate the extinction rate for each species concept.

# For an evolutionary species, the extinction rate is the rate of lineage death.
mu_E = mu

# For a morphospecies, the extinction rate is the sum of the rates of true extinction (mu),
# pseudo-extinction by anagenesis (phi), and pseudo-extinction by bifurcation (0.5*lambda).
mu_M = mu + phi_rate + 0.5 * lambda_rate

# Step 4: Calculate the final ratio and present the result.
# The question asks how greater the extinction rate for a morphospecies is
# compared to that of an evolutionary species. This is the ratio mu_M / mu_E.
multiplicative_factor = mu_M / mu_E

# Print the final equation with all the calculated numbers.
# We explicitly show the values being substituted into the formula for clarity.
print("Calculation Steps:")
print(f"1. Extinction rate for Evolutionary Species (mu_E) = mu = {mu_E}")
print(f"2. Extinction rate for Morphospecies (mu_M) = mu + phi + 0.5*lambda")
print(f"   Substituting derived values (mu={mu}, lambda={lambda_rate}, phi={phi_rate}):")
print(f"   mu_M = {mu} + {phi_rate} + 0.5*{lambda_rate} = {mu_M}")
print("\nFinal Result:")
print(f"The ratio mu_M / mu_E represents the multiplicative factor.")
print(f"Final Equation: mu_M / mu_E = {mu_M} / {mu_E} = {multiplicative_factor}")

# The problem asks for the numeric answer in a specific format.
# We use a file-like object (sys.stdout) to prevent printing the extra characters
# that come with a standard print() of a formatted string, like quotes.
# print(f'<<<{multiplicative_factor}>>>')