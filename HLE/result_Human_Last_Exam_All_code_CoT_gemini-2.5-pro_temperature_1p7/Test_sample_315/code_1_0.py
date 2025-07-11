import math

# Plan:
# 1. Define the physical parameters of the experiment.
# 2. Set representative specific rotation values for D-glucose for different light colors.
# 3. Use the formula for optical rotation (α = [α] * l * c) to calculate the total rotation for each color.
# 4. Print the calculations and results to show the significant difference in rotation angles, which explains the spiral rainbow effect.

# 1. Define experimental parameters
path_length_m = 1.0  # Tube length in meters
# A reasonably high concentration for a visible effect
concentration_g_per_100ml = 50.0

# 2. Define physical constants for D-glucose
# Specific rotation ([α]) is wavelength-dependent. Units are degrees * mL / (g * dm).
# These are approximate values for red, yellow (Na D-line), and blue light.
specific_rotation = {
    "Red": 40.0,
    "Yellow": 52.7,
    "Blue": 100.0,
}

# Convert units for consistency in the formula
# Path length 'l' must be in decimeters (dm)
path_length_dm = path_length_m * 10
# Concentration 'c' must be in g/mL
concentration_g_per_ml = concentration_g_per_100ml / 100.0

print("This problem demonstrates optical rotatory dispersion, where the rotation of polarized light depends on its wavelength (color).")
print("\nCalculating the total rotation for different colors in a D-glucose solution:")
print("------------------------------------------------------------------------")
print(f"Parameters:")
print(f"  Path Length (l): {path_length_m} m = {path_length_dm} dm")
print(f"  Concentration (c): {concentration_g_per_100ml} g/100mL = {concentration_g_per_ml} g/mL")
print("------------------------------------------------------------------------")

# 3. & 4. Calculate and print the rotation for each color
for color, spec_rot in specific_rotation.items():
    # The formula for total rotation angle (alpha) is:
    # alpha = specific_rotation * path_length_dm * concentration_g_per_ml
    total_rotation = spec_rot * path_length_dm * concentration_g_per_ml

    print(f"Calculation for {color.upper()} light:")
    print(f"  Rotation Angle = Specific Rotation * Path Length * Concentration")
    # Outputting each number in the final equation as requested
    print(f"  Angle = {spec_rot} deg*mL/(g*dm) * {path_length_dm} dm * {concentration_g_per_ml} g/mL")
    print(f"  Result: {total_rotation:.1f} degrees\n")

print("The calculations show that blue light rotates significantly more than red light.")
print("This difference causes each color's light to be maximally scattered towards a side observer at a different point along the tube, creating a spiral rainbow effect.")
