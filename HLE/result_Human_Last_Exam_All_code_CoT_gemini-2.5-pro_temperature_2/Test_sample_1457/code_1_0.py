import math

# --- Plan Explanation ---
print("This script determines the lightest element detectable in an EDX setup.")
print("The analysis is based on two main criteria:")
print("1. Excitation: The 30 keV beam must have more energy than the element's K-shell binding energy.")
print("2. Transmission: The element's characteristic X-rays must pass through the 100 µm Beryllium window.")
print("\nWe will calculate the transmission for each element and see which is the lightest element")
print("to pass a reasonable detection threshold (e.g., >1% transmission).\n")


# --- Constants and Setup ---
BEAM_ENERGY_KEV = 30.0
BE_DENSITY_G_CM3 = 1.85  # Density (ρ) of Beryllium
BE_THICKNESS_UM = 100.0
BE_THICKNESS_CM = BE_THICKNESS_UM / 10000.0  # Thickness (t) in cm
TRANSMISSION_THRESHOLD = 0.01  # 1% threshold for detection

# Data for elements: Atomic number (Z), K-alpha X-ray energy (keV),
# K-shell binding energy (keV), and the mass attenuation coefficient (μ)
# of Beryllium at that K-alpha energy (in cm^2/g).
elements_data = {
    'Na': {'Z': 11, 'K_alpha_keV': 1.041, 'K_bind_keV': 1.072, 'mu_Be_cm2_g': 615.0},
    'Mg': {'Z': 12, 'K_alpha_keV': 1.254, 'K_bind_keV': 1.303, 'mu_Be_cm2_g': 357.0},
    'Si': {'Z': 14, 'K_alpha_keV': 1.740, 'K_bind_keV': 1.839, 'mu_Be_cm2_g': 158.4},
    'Ca': {'Z': 20, 'K_alpha_keV': 3.692, 'K_bind_keV': 4.038, 'mu_Be_cm2_g': 27.2}
}

print(f"--- Analysis Parameters ---")
print(f"Beam Energy: {BEAM_ENERGY_KEV} keV")
print(f"Beryllium Window Thickness (t): {BE_THICKNESS_CM} cm")
print(f"Beryllium Density (ρ): {BE_DENSITY_G_CM3} g/cm³")
print(f"Detection Threshold: {TRANSMISSION_THRESHOLD * 100}% transmission")
print("---------------------------\n")

detectable_elements = []

# Sort elements by atomic number to evaluate the lightest ones first
sorted_elements = sorted(elements_data.items(), key=lambda item: item[1]['Z'])

for name, data in sorted_elements:
    print(f"--- Analyzing: {name} (Z={data['Z']}) ---")

    # 1. Check Excitation
    is_excitable = BEAM_ENERGY_KEV > data['K_bind_keV']
    print(f"1. Excitation Check:")
    print(f"   Is Beam Energy ({BEAM_ENERGY_KEV} keV) > Binding Energy ({data['K_bind_keV']} keV)? {'Yes.' if is_excitable else 'No.'}")

    if not is_excitable:
        print("   Element cannot be excited. Not detectable.\n")
        continue

    # 2. Calculate Transmission
    mu = data['mu_Be_cm2_g']
    rho = BE_DENSITY_G_CM3
    t = BE_THICKNESS_CM
    
    # Calculate the exponent for the Beer-Lambert law
    attenuation_factor = mu * rho * t
    transmission = math.exp(-attenuation_factor)

    print("2. Transmission Check:")
    print(f"   X-ray (K-alpha) Energy: {data['K_alpha_keV']} keV")
    print(f"   The transmission formula is T = exp(-μ * ρ * t)")
    # Print the equation with all the numbers
    print(f"   T = exp(-{mu:.1f} * {rho:.2f} * {t:.4f})")
    print(f"   T = exp(-{attenuation_factor:.4f}) = {transmission:.6f}")
    print(f"   Transmission Percentage: {transmission * 100:.4f}%")
    
    is_detectable = transmission > TRANSMISSION_THRESHOLD
    print(f"   Is Transmission > {TRANSMISSION_THRESHOLD * 100}%? {'Yes, this element is likely detectable.' if is_detectable else 'No, this element is likely not detectable.'}")

    if is_detectable:
        detectable_elements.append(name)
    
    print("-" * 28 + "\n")

# --- Conclusion ---
print("--- Final Conclusion ---")
if detectable_elements:
    lightest_detectable_element = detectable_elements[0]
    print(f"The elements with calculated transmission > 1% are: {', '.join(detectable_elements)}.")
    print(f"The lightest element among these is {lightest_detectable_element} ({elements_data[lightest_detectable_element]['Z']}).")
    if lightest_detectable_element == 'Si':
        print("\nAnswer: Based on the analysis, Silicon (Si) is the lightest element from the list that can be reliably detected.")
else:
    print("\nBased on a 1% transmission threshold, none of the light elements would be detectable with this setup.")
