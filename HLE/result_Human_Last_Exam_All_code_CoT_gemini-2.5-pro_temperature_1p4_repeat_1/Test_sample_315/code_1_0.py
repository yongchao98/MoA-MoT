import numpy as np

# --- Setup Parameters ---
# D-glucose solution properties
# Concentration in g/mL. A high concentration is used for a dramatic effect.
concentration_g_ml = 0.7 
# Specific rotation for Sodium D-line (589 nm) in degrees*mL/(g*dm)
specific_rotation_Na = 52.7
# Wavelength of Sodium D-line in nm
lambda_Na = 589.0

# Tube dimensions
# Length in meters, converted to decimeters (dm) for the formula
length_m = 1.0
length_dm = length_m * 10

# --- Physics Calculation ---
# This phenomenon is called Optical Rotatory Dispersion (ORD).
# The specific rotation [a] is approximately proportional to 1/lambda^2.
# We can find the proportionality constant k from the known value for the Sodium D-line.
# [a] = k / lambda^2  =>  k = [a]_Na * lambda_Na^2
k = specific_rotation_Na * (lambda_Na**2)

# Wavelengths for different colors in nm
colors = {
    "Violet": 410,
    "Blue": 470,
    "Green": 520,
    "Yellow": 580,
    "Orange": 610,
    "Red": 680,
}

print(f"Calculating the total rotation for different colors of light in a {length_m}m tube of D-glucose solution.")
print(f"Concentration = {concentration_g_ml} g/mL, Length = {length_dm} dm\n")
print("-" * 55)
print(f"{'Color':<10} | {'Wavelength (nm)':<18} | {'Total Rotation (degrees)':<22}")
print("-" * 55)

# Calculate and print rotation for each color
for color, wavelength in colors.items():
    # Calculate specific rotation for this wavelength based on ORD
    specific_rotation_lambda = k / (wavelength**2)
    
    # Calculate total rotation angle using the formula: angle = [a] * L * c
    total_rotation = specific_rotation_lambda * length_dm * concentration_g_ml
    
    # Calculate number of full 360-degree turns
    turns = total_rotation / 360.0
    
    print(f"{color:<10} | {wavelength:<18} | {total_rotation:<22.0f} ({turns:.1f} turns)")

print("-" * 55)
print("\nConclusion:")
print("As the calculations show, different colors of light rotate by vastly different amounts.")
print("Violet light completes more than twice as many turns as red light over the same distance.")
print("This differential rotation, combined with light scattering, creates a visible rainbow effect that spirals down the tube.")