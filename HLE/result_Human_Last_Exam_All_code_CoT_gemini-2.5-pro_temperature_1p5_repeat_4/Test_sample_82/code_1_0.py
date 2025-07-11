import sys
# Redirect stdout to a variable to prevent printing the final answer marker
from io import StringIO
original_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

# --- Problem Setup ---
# For simplicity, let's assume the total Phenotypic Variance (Vp) is 1.0.
Vp = 1.0
# We are given the Broad-Sense Heritability (H^2).
H2 = 0.5

# --- Calculations ---
# From the definition H^2 = Vg / Vp, we can find the total Genetic Variance (Vg).
Vg = H2 * Vp

print("--- Step-by-Step Analysis ---")
print(f"Given: Total Phenotypic Variance (Vp) = {Vp}")
print(f"Given: Broad-Sense Heritability (H^2) = {H2}")
print(f"Calculated: Total Genetic Variance (Vg = H^2 * Vp) = {H2} * {Vp} = {Vg}\n")

print("--- Analysis of Statement A ---")
print("A Polygenic Score (PGS) is a predictor based on genetic data.")
print("Its predictive power is fundamentally limited by the total contribution of genes to the phenotype, which is Vg.")
print(f"The maximum variance a PGS can explain is {Vg}, or {Vg/Vp*100}%.")
print("Statement A says the PGS cannot explain MORE than 50%. This is TRUE by definition.\n")

print("--- Analysis of Statements B and C ---")
print("Genetic Variance (Vg) is composed of Additive (Va), Dominance (Vd), and Epistatic (Vi) components.")
print(f"Equation: Vg = Va + Vd + Vi = {Vg}")
print("A standard PGS from GWAS sums linear effects and best captures Additive Variance (Va).")
print("The theoretical maximum variance a PGS can explain is Narrow-Sense Heritability (h^2 = Va / Vp).")

print("\nLet's model a realistic scenario for a complex trait:")
# Assume some variance is non-additive (dominance, epistasis)
Va = 0.35
Vd = 0.1
Vi = 0.05
calculated_Vg = Va + Vd + Vi
print(f"  Let Va = {Va}")
print(f"  Let Vd = {Vd}")
print(f"  Let Vi = {Vi}")
print(f"  Check: Va + Vd + Vi = {Va} + {Vd} + {Vi} = {calculated_Vg}. This equals our Vg.")

# Calculate the resulting narrow-sense heritability
h2 = Va / Vp
print(f"\nIn this realistic case, h^2 = Va / Vp = {Va} / {Vp} = {h2}")
print(f"The PGS would approach an explained variance of {h2*100}%, NOT {H2*100}%.")
print("Therefore, Statement B (will approach 50%) is FALSE.")
print("And Statement C (will NOT approach 50% due to non-linear effects) is TRUE.\n")

print("--- Analysis of Statement D ---")
print("Epigenetic effects are generally considered part of Environmental Variance (Ve).")
print("The existence of epigenetic effects does not dictate the composition of Vg.")
print("We could have a scenario where Vd=0 and Vi=0, making Vg = Va.")
print(f"  E.g., Va = {Vg}, Vd = 0, Vi = 0. In this case, h^2 = Va / Vp = {Vg} / {Vp} = {H2}")
print(f"Even with epigenetics, h^2 could equal {H2}. So epigenetics does not force h^2 to be less than 0.5.")
print("Therefore, Statement D is FALSE.\n")

print("--- FINAL CONCLUSION ---")
print("Statements A and C are the only ones that are necessarily true.")
print("This corresponds to answer choice E.")

# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())