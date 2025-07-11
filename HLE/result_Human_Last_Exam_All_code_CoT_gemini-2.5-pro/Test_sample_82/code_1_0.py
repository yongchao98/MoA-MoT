import numpy as np

# Set a random seed for reproducible results
np.random.seed(42)

# --- 1. Define Population and Genetic Parameters ---
# Use a large number of individuals for stable variance estimates
n_individuals = 50000
n_snps = 100

# Heritability parameters as given in the problem
H2 = 0.5  # Broad-sense heritability (Total Genetic Variance / Total Phenotypic Variance)

# To test the statements, especially C, we must assume some non-additive genetic variance exists.
# This means narrow-sense heritability (h²) must be less than broad-sense heritability (H²).
# Let's set h² to a value less than H².
h2 = 0.4  # Narrow-sense heritability (Additive Genetic Variance / Total Phenotypic Variance)

# The remaining portion of genetic variance is non-additive (in this case, epistatic)
epistatic_variance_prop = H2 - h2  # This will be 0.1

# The rest of the phenotypic variance is due to the environment
environmental_variance_prop = 1 - H2  # This will be 0.5

# --- 2. Simulate Genotypes for the Population ---
# Assume SNPs are independent and simulate genotypes (0, 1, or 2 minor alleles)
mafs = np.random.uniform(0.05, 0.5, n_snps)
genotypes = np.random.binomial(2, mafs, size=(n_individuals, n_snps))

# --- 3. Simulate Genetic and Environmental Components of the Phenotype ---

# a) Additive Genetic Component (A)
# A standard PGS is a linear model based on summing the effects of individual SNPs.
additive_effects = np.random.randn(n_snps)
raw_additive_values = genotypes @ additive_effects
# Scale the additive component to have the desired variance (h2)
scaled_additive_values = raw_additive_values / np.std(raw_additive_values) * np.sqrt(h2)

# b) Non-additive (Epistatic) Genetic Component (I)
# We model a simple gene-gene interaction between the first two SNPs.
interaction_term = (genotypes[:, 0] == 2) & (genotypes[:, 1] == 2)
# Scale the epistatic component to have the desired variance (H2 - h2)
scaled_epistatic_values = interaction_term / np.std(interaction_term) * np.sqrt(epistatic_variance_prop)

# c) Total Genetic Value (G = A + I)
total_genetic_value = scaled_additive_values + scaled_epistatic_values

# d) Environmental Component (E)
environmental_values = np.random.randn(n_individuals) * np.sqrt(environmental_variance_prop)

# --- 4. Simulate the Final Phenotype (P = G + E) ---
phenotype = total_genetic_value + environmental_values
# Center the phenotype to have a mean of 0 for easier variance calculation
phenotype = phenotype - np.mean(phenotype)

# --- 5. Analyze Results and Evaluate Statements ---
var_p = np.var(phenotype)
var_g = np.var(total_genetic_value)
var_a = np.var(scaled_additive_values)

# The variance explained by a perfect, ideal PGS is the squared correlation
# between the true additive genetic values and the final phenotype.
r2_perfect_pgs = (np.corrcoef(scaled_additive_values, phenotype)[0, 1])**2

print("--- Simulation Results ---")
print(f"Target Broad-Sense Heritability (H²): {H2}")
print(f"Target Narrow-Sense Heritability (h²): {h2}")
print("-" * 30)
print(f"Simulated Total Genetic Variance (Vg/Vp): {var_g/var_p:.4f}")
print(f"Simulated Additive Genetic Variance (Va/Vp): {var_a/var_p:.4f}")
print("-" * 30)

print("\n--- Statement Analysis ---")
print("\nStatement A: The polygenic score can not explain more than 50% of the variance.")
print("The maximum variance explainable by genetics is the total genetic variance (Vg/Vp).")
print(f"In our simulation, Vg/Vp is {var_g/var_p:.4f}, which is equal to H² (0.5).")
print("A PGS, being based on genetics, cannot exceed this limit. Statement A is TRUE.")

print("\nStatement C: A linear PGS will not approach a variance explained of 50% due to non-linear effects.")
print("A perfect linear PGS can only capture the additive variance (Va/Vp).")
print(f"The variance explained by our simulated ideal PGS is R² = {r2_perfect_pgs:.4f}.")
print(f"This value is approximately h² ({h2}), not H² ({H2}), because non-additive (epistatic) effects were included.")
print("The gap is due to these non-linear effects. Statement C is TRUE.")

print("\n--- Final Conclusion ---")
print("Both statements A and C are necessarily true.")
print("The correct option is E, which states that only choices A and C are correct.")
<<<E>>>