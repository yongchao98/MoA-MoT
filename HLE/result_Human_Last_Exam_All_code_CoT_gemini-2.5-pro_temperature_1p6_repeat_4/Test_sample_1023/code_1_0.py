# Step 1: Define the extinction rate for an evolutionary species (μ_e).
# This rate is determined by the true biological extinction of a lineage.
# Let's denote the rate of true biological extinction as the variable 'μ'.
mu_e = "μ"
print(f"The extinction rate for an evolutionary species, μ_e, is the rate of true extinction: {mu_e}")

# Step 2: Define the extinction rate for a morphospecies (μ_m).
# This rate is the sum of all processes that lead to the disappearance of a morphospecies.
# These processes are:
#  - True biological extinction (rate μ).
#  - Pseudo-extinction from anagenesis (rate α).
#  - Pseudo-extinction from bifurcating speciation (rate 0.5 * λ, where λ is the branching rate).
# The total rate is μ_m = μ + α + 0.5 * λ.
# We can group the artifactual extinctions into a single term: μ_pseudo = α + 0.5 * λ.
mu_pseudo_formula = "α + 0.5 * λ"
print(f"\nThe total rate of extinction for a morphospecies, μ_m, is the sum of true extinction (μ) and pseudo-extinction (μ_pseudo).")
print(f"The rate of pseudo-extinction is μ_pseudo = {mu_pseudo_formula}.")
print("So, the total extinction rate is μ_m = μ + μ_pseudo.")

# Step 3: Apply the central assumption from the problem statement.
# The assumption is that "all the processes that affect them occur at the same rates".
# This is interpreted to mean that the rate of true biological extinction is equal to the rate of pseudo-extinction.
print("\nThe key assumption is that the rate of true extinction is equal to the rate of pseudo-extinction.")
print("Assumption: μ = μ_pseudo")

# Step 4: Calculate the final extinction rate for a morphospecies based on the assumption.
# By substituting μ_pseudo with μ, we get:
mu_m_final_formula = "μ + μ"
print(f"\nUsing this assumption, the morphospecies extinction rate becomes: μ_m = {mu_m_final_formula} = 2 * μ.")

# Step 5: Calculate the final ratio to answer the question.
# The question asks how much greater the extinction rate for a morphospecies is compared to an evolutionary species.
# This is the ratio μ_m / μ_e.
print("\nWe now compute the ratio of the two extinction rates:")
numerator = "2 * μ"
denominator = "μ"
result = 2
print(f"Ratio = μ_m / μ_e = ({numerator}) / ({denominator}) = {result}")
print("\nThis means the extinction rate for a morphospecies is twice the extinction rate for an evolutionary species.")
print("\nFinal equation with numbers shown:")
print(f"2*μ / 1*μ = 2")

<<<2>>>