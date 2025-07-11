import xraylib
import math

# Plan:
# 1. Identify the physics: To be detected, an element's characteristic X-rays must pass through the detector's 100 µm Beryllium window.
# 2. Formulate the problem: We must calculate the transmission percentage for the K-alpha X-rays of each candidate light element (Na, Mg, Si, Ca).
# 3. The lightest element with a transmission percentage high enough to be measured (e.g., >1%) is the answer.
# 4. Use the Beer-Lambert law: Transmission = exp(-μ * ρ * x), where μ is the mass attenuation coefficient, ρ is the density, and x is the thickness.
# 5. Use the 'xraylib' library for accurate physical data (X-ray energies, attenuation coefficients, density).

# --- Data and Constants ---
# Window parameters
WINDOW_MATERIAL_SYMBOL = 'Be'
WINDOW_MATERIAL_Z = xraylib.SymbolToAtomicNumber(WINDOW_MATERIAL_SYMBOL)
WINDOW_THICKNESS_UM = 100
WINDOW_THICKNESS_CM = WINDOW_THICKNESS_UM * 1e-4  # Convert µm to cm
WINDOW_DENSITY_G_CM3 = xraylib.ElementDensity(WINDOW_MATERIAL_Z)

# Elements to analyze from the answer choices (name: atomic number)
# Sorted by atomic number to find the "lightest" first.
ELEMENTS_TO_CHECK = {
    'Na': 11,
    'Mg': 12,
    'Si': 14,
    'Ca': 20
}

# A reasonable threshold for a signal to be considered "detectable"
DETECTION_THRESHOLD_PERCENT = 1.0

# --- Calculation and Output ---
print("Analyzing which is the lightest element detectable with a {} µm {} window.\n".format(WINDOW_THICKNESS_UM, WINDOW_MATERIAL_SYMBOL))

# --- Example Calculation Showing Each Number ---
# As requested, showing the full equation for one element (Silicon)
si_z = ELEMENTS_TO_CHECK['Si']
si_energy_keV = xraylib.LineEnergy(si_z, xraylib.KA_LINE)
si_mu_cm2_g = xraylib.CS_Total(WINDOW_MATERIAL_Z, si_energy_keV)
si_transmission_percent = 100 * math.exp(-si_mu_cm2_g * WINDOW_DENSITY_G_CM3 * WINDOW_THICKNESS_CM)

print("To find the answer, we use the Beer-Lambert law: T(%) = 100 * exp(-μ * ρ * x)")
print("\nExample calculation for Silicon (Si):")
print("Here, we plug in each number for the equation:")
print("μ (mass attenuation coefficient of Be for Si Kα X-rays) = {:.1f} cm^2/g".format(si_mu_cm2_g))
print("ρ (density of Be) = {:.3f} g/cm^3".format(WINDOW_DENSITY_G_CM3))
print("x (thickness of Be window) = {} cm".format(WINDOW_THICKNESS_CM))
print("\nFinal equation for Si:")
print("T(%) = 100 * exp(-{:.1f} * {:.3f} * {})".format(si_mu_cm2_g, WINDOW_DENSITY_G_CM3, WINDOW_THICKNESS_CM))
print("T(%) for Si = {:.4f}%\n".format(si_transmission_percent))


# --- Results for all elements ---
print("-" * 80)
print("{:<10} | {:<4} | {:<12} | {:<12} | {:<18} | {:<10}".format('Element', 'Z', 'Energy (keV)', 'μ (cm^2/g)', 'Transmission (%)', 'Detectable?'))
print("-" * 80)

lightest_element_found = None
final_answer_choice = None

for name, z in sorted(ELEMENTS_TO_CHECK.items(), key=lambda item: item[1]):
    # Get the K-alpha X-ray energy for the element in keV
    ka_energy_keV = xraylib.LineEnergy(z, xraylib.KA_LINE)

    # Get the mass attenuation coefficient (μ) of Beryllium for that energy
    mu_cm2_g = xraylib.CS_Total(WINDOW_MATERIAL_Z, ka_energy_keV)

    # Calculate transmission percentage
    exponent = -mu_cm2_g * WINDOW_DENSITY_G_CM3 * WINDOW_THICKNESS_CM
    transmission_percent = 100 * math.exp(exponent)
    
    # Check if it meets the detection threshold
    is_detected = "Yes" if transmission_percent > DETECTION_THRESHOLD_PERCENT else "No"
    
    print("{:<10} | {:<4} | {:<12.3f} | {:<12.1f} | {:<18.4f} | {:<10}".format(name, z, ka_energy_keV, mu_cm2_g, transmission_percent, is_detected))
    
    # Store the first element that is detected
    if is_detected == "Yes" and lightest_element_found is None:
        lightest_element_found = name
        # Map element name back to the provided answer choices
        answer_map = {'Na': 'B', 'Mg': 'A', 'Si': 'E', 'Ca': 'D'}
        final_answer_choice = answer_map[name]

print("-" * 80)

if lightest_element_found:
    print("\nConclusion: Based on a {:.1f}% detection threshold, the lightest detectable element is {}.".format(DETECTION_THRESHOLD_PERCENT, lightest_element_found))
    print("This corresponds to answer choice {}.".format(final_answer_choice))
else:
    print("\nConclusion: None of the listed elements have a transmission greater than {:.1f}%.".format(DETECTION_THRESHOLD_PERCENT))
