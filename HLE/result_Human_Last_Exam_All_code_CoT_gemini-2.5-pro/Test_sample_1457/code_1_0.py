import xraylib
import math

# Plan:
# 1. The main obstacle for detecting light elements in EDX is the absorption of their low-energy X-rays by the detector's window.
# 2. We will calculate the transmission of the characteristic Kα X-rays for each candidate element (Na, Mg, Si, Ca) through the specified Be window.
# 3. The calculation uses the Beer-Lambert Law: Transmission = exp(-μ * ρ * x), where:
#    - μ is the mass attenuation coefficient of Be for the specific X-ray energy.
#    - ρ is the density of Be.
#    - x is the thickness of the Be window.
# 4. We will use the 'xraylib' library for accurate physical data (energies, μ, ρ).
# 5. The lightest element with a transmission percentage high enough to be detected (e.g., > 1%) is the correct answer.

# --- Parameters and Constants ---
WINDOW_MATERIAL = 'Be'
WINDOW_THICKNESS_UM = 100
WINDOW_THICKNESS_CM = WINDOW_THICKNESS_UM * 1e-4  # Convert µm to cm

# Candidate elements from the answer choices and their atomic numbers (Z)
CANDIDATE_ELEMENTS = {'Na': 11, 'Mg': 12, 'Si': 14, 'Ca': 20}

# Get the density of the Beryllium window using xraylib
# Z of Beryllium is 4
WINDOW_DENSITY_G_CM3 = xraylib.ElementDensity(4)

print("Analyzing which lightest element is detectable with a 100 µm Beryllium window.\n")

# Loop through elements, sorted by atomic number (lightest to heaviest)
for name, Z in sorted(CANDIDATE_ELEMENTS.items(), key=lambda item: item[1]):
    print(f"--- Checking Element: {name} (Z={Z}) ---")

    # Get the energy of the K-alpha characteristic X-ray line in keV
    energy_keV = xraylib.LineEnergy(Z, xraylib.KA_LINE)
    print(f"Characteristic Kα X-ray Energy: {energy_keV:.4f} keV")

    # Get the mass attenuation coefficient (μ) of Be for this energy
    mu_cm2_g = xraylib.CS_Total_CP(WINDOW_MATERIAL, energy_keV)
    print(f"Mass Attenuation Coeff. (μ) of Be at this energy: {mu_cm2_g:.2f} cm²/g")

    # Calculate the transmission using the Beer-Lambert equation
    exponent = -mu_cm2_g * WINDOW_DENSITY_G_CM3 * WINDOW_THICKNESS_CM
    transmission = math.exp(exponent)
    transmission_percent = transmission * 100

    # Print the full calculation for clarity
    print(f"Transmission = exp(-μ * ρ * x)")
    print(f"             = exp(-{mu_cm2_g:.2f} * {WINDOW_DENSITY_G_CM3:.3f} * {WINDOW_THICKNESS_CM})")
    print(f"             = exp({exponent:.4f})")
    print(f"Resulting Transmission: {transmission_percent:.6f} %")

    if transmission_percent > 1.0:
        print("Conclusion: Detectable (transmission > 1%).\n")
    else:
        print("Conclusion: Effectively undetectable (transmission < 1%).\n")

print("Based on the calculations, the lightest element from the list with a significant transmission (>1%) is Calcium (Ca).")
<<<D>>>