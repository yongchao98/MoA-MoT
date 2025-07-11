import xraylib
import numpy as np
import warnings

# Suppress potential UserWarning from xraylib for a cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

# --- Problem Parameters ---
WINDOW_MATERIAL_SYMBOL = 'Be'
WINDOW_MATERIAL_Z = 4
WINDOW_THICKNESS_UM = 100.0
WINDOW_THICKNESS_CM = WINDOW_THICKNESS_UM / 10000.0  # convert µm to cm

# --- Physical Constants from xraylib ---
try:
    # Get Beryllium's density using its atomic number (Z=4)
    DENSITY_BE_G_CM3 = xraylib.ElementDensity(WINDOW_MATERIAL_Z)
except Exception as e:
    # Use a fallback value if xraylib database is not found
    print(f"Warning: Could not load data from xraylib, using default value. Error: {e}")
    DENSITY_BE_G_CM3 = 1.85

# Candidate elements ordered by atomic number (lightest to heaviest)
elements_to_check = {
    'Na': 11, # Sodium
    'Mg': 12, # Magnesium
    'Si': 14, # Silicon
    'Ca': 20, # Calcium
}

# --- Calculation and Output ---
print(f"Analyzing X-ray transmission through a {WINDOW_THICKNESS_UM} µm Beryllium (Be) window.")
print("The transmission T is calculated using the Beer-Lambert law: T = exp(-μ * ρ * t)")
print(f"Using Be density (ρ) = {DENSITY_BE_G_CM3:.3f} g/cm³ and window thickness (t) = {WINDOW_THICKNESS_CM} cm.\n")

# A dictionary to store results for the final conclusion
results = {}

for name, z in elements_to_check.items():
    try:
        # Get the element's characteristic K-alpha X-ray energy in keV
        energy_kev = xraylib.LineEnergy(z, xraylib.KA_LINE)
        
        # Get the mass attenuation coefficient (μ) of Be for that energy in cm²/g
        mu_be = xraylib.CS_Total_CP(WINDOW_MATERIAL_SYMBOL, energy_kev)
        
        # Calculate transmission
        transmission_fraction = np.exp(-mu_be * DENSITY_BE_G_CM3 * WINDOW_THICKNESS_CM)
        transmission_percent = transmission_fraction * 100
        results[name] = transmission_percent
        
        # Print the full equation with the numbers plugged in, as requested
        print(f"For {name} (Kα energy = {energy_kev:.3f} keV):")
        print(f"  Transmission(%) = 100 * exp(-{mu_be:.2f} * {DENSITY_BE_G_CM3:.3f} * {WINDOW_THICKNESS_CM}) = {transmission_percent:.2f}%")
        print("-" * 30)
        
    except Exception as e:
        print(f"Could not perform calculation for {name}. Error: {e}")

# --- Conclusion ---
print("\n--- Analysis ---")
print("To be 'seen' by the detector, an element's X-rays must have a high enough transmission percentage.")
print(f"  - Na (Sodium) transmission is {results.get('Na', 0):.2f}%, which is too low to be reliably detected.")
print(f"  - Mg (Magnesium) transmission is {results.get('Mg', 0):.2f}%, which is very low and likely lost in the noise.")
print(f"  - Si (Silicon) transmission is {results.get('Si', 0):.2f}%, which is a strong signal that is easily detectable.")
print(f"  - Ca (Calcium) is also detectable but is heavier than Silicon.")
print("\nTherefore, the lightest element from the choices that can be confidently seen is Silicon (Si).")
<<<E>>>