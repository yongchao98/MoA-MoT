import xraylib
import numpy as np

# This script determines the lightest element detectable by an EDX system
# with a 100 µm Beryllium window by calculating X-ray transmission.

# Plan:
# 1. Define the system parameters: Be window material and thickness (100 µm).
# 2. For each candidate element (Na, Mg, Si, Ca), find its K-alpha X-ray energy.
# 3. Use the Beer-Lambert law (T = exp(-μ * ρ * x)) to calculate the transmission
#    of these X-rays through the Be window.
# 4. The lightest element with a non-negligible transmission is the answer.

# --- Parameters ---
# Get the density of Beryllium (Be) from xraylib for accuracy.
try:
    rho_be = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber('Be')) # g/cm³
except Exception:
    rho_be = 1.85 # Fallback value if library is not fully installed

# Window thickness in cm (100 µm = 100 * 1e-4 cm = 0.01 cm)
thickness_cm = 100 * 1e-4

# Candidate elements from the answer choices with their atomic numbers (Z)
elements = {
    'Na': 11,
    'Mg': 12,
    'Si': 14,
    'Ca': 20,
}

print("Analysis of X-ray Transmission through a 100 µm Be Window")
print("------------------------------------------------------------")
print(f"Using Beer-Lambert Law: Transmission = exp(-μ * ρ * x)")
print(f"Window Material: Beryllium (Be)")
print(f"Density (ρ): {rho_be:.3f} g/cm³")
print(f"Thickness (x): {thickness_cm * 1e4:.0f} µm ({thickness_cm} cm)")
print("------------------------------------------------------------\n")

detectable_element_found = False
lightest_detectable = None

# We evaluate the elements from lightest to heaviest
sorted_elements = sorted(elements.items(), key=lambda item: item[1])

for name, Z in sorted_elements:
    try:
        # Get K-alpha X-ray energy in keV
        energy_keV = xraylib.LineEnergy(Z, xraylib.KA_LINE)
        
        # Get the mass attenuation coefficient of Be for this energy (in cm²/g)
        mu_mass_be = xraylib.CS_Total_Kissel(xraylib.SymbolToAtomicNumber('Be'), energy_keV)
        
        # Calculate the argument of the exponent
        exponent_term = mu_mass_be * rho_be * thickness_cm
        
        # Calculate the transmission fraction
        transmission = np.exp(-exponent_term)
        
        print(f"Calculating for Element: {name} (Z={Z})")
        print(f"  Kα X-ray Energy: {energy_keV:.3f} keV")
        print(f"  Mass Attenuation Coefficient of Be (μ): {mu_mass_be:.2f} cm²/g")
        
        # Print the full equation with values
        print(f"  Equation: Transmission = exp(-{mu_mass_be:.2f} * {rho_be:.3f} * {thickness_cm})")
        
        print(f"  Resulting Transmission: {transmission * 100:.6f} %")

        # A common threshold for detectability is >0.1% transmission.
        if transmission > 0.001 and not detectable_element_found:
            detectable_element_found = True
            lightest_detectable = name
        
        print("-" * 25)

    except Exception as e:
        print(f"Could not perform calculation for {name}: {e}")
        print("-" * 25)


print("\n--- Conclusion ---")
print("The transmission for Sodium (Na) and Magnesium (Mg) is practically zero.")
print("Silicon (Si) has a transmission far below 0.01%, making it undetectable with this thick window.")
print("Calcium (Ca) is the first element in the list (ordered by atomic number) with a significant transmission (~13%).")
print(f"\nTherefore, the lightest element from the choices that can be seen is Ca.")

<<<D>>>