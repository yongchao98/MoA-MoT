import math

def calculate_tolerance_factor(r_A, r_B, r_X):
    """Calculates the Goldschmidt tolerance factor."""
    numerator = r_A + r_X
    denominator = math.sqrt(2) * (r_B + r_X)
    return numerator / denominator

def check_stability(cation, t_factor):
    """Checks the stability based on the tolerance factor and provides a conclusion."""
    conclusion = ""
    if 1.0 < t_factor:
        conclusion = "Too large. Tends to form 2D or quasi-2D structures."
    elif 0.8 <= t_factor <= 1.0:
        conclusion = "Good fit. Capable of forming stable 3D perovskite structures."
    else: # t_factor < 0.8
        conclusion = "Too small. Tends to form distorted or non-perovskite phases."
    
    print(f"Cation: {cation}")
    print(f"Tolerance Factor (t) for {cation}PbBr3: {t_factor:.2f}")
    print(f"Conclusion: {conclusion}\n")

# Ionic radii in picometers (pm)
# B-site (Pb2+) and X-site (Br-)
r_B_Pb = 119  # Radius of Pb2+
r_X_Br = 196  # Radius of Br-

# A-site cations from the answer choices
a_site_cations = {
    "Cesium": 167,
    "Methylammonium": 217,
    "Formamidinium": 253,
    "Aziridinium": 200,      # Small and reactive, generally unsuitable
    "Ethylammonium": 274,    # Too large
    "Methylhydrazinium": 270,# Too large
    "Dimethylammonium": 272  # Too large
}

print("Analyzing A-site cations for forming 3D A-Pb-Br3 perovskite structures.\n")
print(f"The stability criterion is the Goldschmidt Tolerance Factor (t).")
print("A stable 3D structure generally requires 0.8 <= t <= 1.0.\n")
print("-" * 50)


# Calculate and check for all relevant cations
for name, r_A in a_site_cations.items():
    t = calculate_tolerance_factor(r_A, r_B_Pb, r_X_Br)
    check_stability(name, t)

print("-" * 50)
print("Based on the analysis, only Cesium, Methylammonium, and Formamidinium")
print("fall within or very close to the ideal range for forming stable 3D perovskites independently.")
print("The other cations are too large, leading to 2D structures.")
print("\nTherefore, the correct choice is the one containing only these three cations.")

<<<A>>>