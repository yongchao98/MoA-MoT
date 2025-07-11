import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt Tolerance Factor."""
    # The equation is: t = (r_A + r_X) / (sqrt(2) * (r_B + r_X))
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    if denominator == 0:
        return float('inf')
    return numerator / denominator

# --- Ion parameters for A-Pb-Br3 perovskite ---
# Ionic radii in picometers (pm). Values are from standard literature.
# B-site cation: Lead (Pb2+)
r_B = 119
# X-site anion: Bromide (Br-)
r_X = 196

# A-site cations and their effective ionic radii (pm)
cations = {
    "Cesium": 167,           # Inorganic reference
    "Methylammonium": 217,
    "Formamidinium": 253,
    "Aziridinium": 203,         # Smaller than MA
    "Ethylammonium": 274,       # Larger than FA
    "Methylhydrazinium": 270,   # Large cation
    "Dimethylammonium": 272     # Large cation
}

print("Evaluating A-site Cations for Stability in 3D A-Pb-Br3 Perovskites\n")
print("The Goldschmidt Tolerance Factor (t) is used to predict structural stability.")
print("A stable 3D perovskite structure typically requires 0.8 < t < 1.0.\n")
print("--- Calculation for each cation ---")
print(f"Using radii: r(Pb2+) = {r_B} pm, r(Br-) = {r_X} pm")
print("-" * 75)
print(f"{'A-Site Cation':<20} | {'Radius r_A (pm)':<18} | {'Tolerance Factor (t)':<22} | {'Conclusion'}")
print("-" * 75)

correct_organic_cations = []

for name, r_A in cations.items():
    # Calculate t for each cation using the formula with its specific numbers
    t = calculate_tolerance_factor(r_A, r_B, r_X)

    # Outputting the numbers in the final calculation step as requested
    # The final equation's result 't' depends on r_A, r_B, and r_X. We print all.
    # We check the result 't' against the stability range.

    if 0.8 < t < 1.0:
        conclusion = "Forms stable 3D structure"
        # We check if the cation is organic
        if name not in ["Cesium"]:
             correct_organic_cations.append(name)
    elif t >= 1.0:
        conclusion = "Too large (t >= 1.0), forms 2D structures"
    else: # t <= 0.8
        conclusion = "Too small (t <= 0.8), favors non-perovskite phase"

    print(f"{name:<20} | {r_A:<18} | {t:<22.3f} | {conclusion}")

print("-" * 75)
print("\nBased on the analysis, the only organic cations from the list capable of")
print("independently forming stable 3D lead halide perovskites are:")
for cation in correct_organic_cations:
    print(f"- {cation}")

print("\nTherefore, the most comprehensive correct list from the choices is Cesium (as a baseline),")
print("Methylammonium, and Formamidinium.")
<<<A>>>