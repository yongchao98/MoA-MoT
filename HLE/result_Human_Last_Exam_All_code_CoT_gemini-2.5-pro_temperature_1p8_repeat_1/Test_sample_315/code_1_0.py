import math

# Plan: To determine the appearance of the tube, we need to understand
# how D-glucose solution interacts with light. D-glucose is optically active,
# meaning it rotates the plane of polarized light. This effect, known as
# optical rotation, is wavelength-dependent (Optical Rotatory Dispersion).
# Shorter wavelengths (blue) are rotated more than longer wavelengths (red).
#
# As light travels down the tube, the plane of polarization for each color
# rotates at a different rate. This creates a helical "fan" of polarization
# planes. Light scatters off the solution, and the intensity seen from the side
# depends on the orientation of the polarization plane. Because of the different
# rotation rates, different colors will be most intensely scattered towards the
# observer at different points along the tube. This results in a spiraling
# rainbow pattern.
#
# The code below will calculate the total rotation for red, green, and blue
# light to demonstrate this dispersion effect numerically.

# --- Parameters ---
# Length of the tube in decimeters (1 meter = 10 dm)
length_dm = 10.0
# A plausible concentration for D-glucose solution in g/mL (e.g., 20% w/v)
concentration_g_ml = 0.20

# Approximate Cauchy's equation constants for D-glucose.
# The specific rotation [α] is given by: [α] = A + B / λ^2
# A is a constant in deg*mL/(g*dm)
A_const = 8.35
# B is a constant in deg*mL*nm^2/(g*dm)
B_const = 1.48e6

# Wavelengths (λ) for representative colors in nanometers (nm)
wavelengths_nm = {
    "Red": 650.0,
    "Green": 550.0,
    "Blue": 450.0
}

print("This script calculates the total optical rotation for different colors of light")
print("passing through a D-glucose solution, demonstrating optical rotatory dispersion.")
print("\nThe formula for total rotation is:")
print("Total Rotation = (A + B / Wavelength^2) * Length * Concentration\n")
print(f"Using: Length = {length_dm} dm, Concentration = {concentration_g_ml} g/mL\n")
print("-" * 70)

# Calculate and display the final equation for each color
for color, wavelength in wavelengths_nm.items():
    # Specific rotation [α] in deg*mL/(g*dm)
    specific_rotation = A_const + (B_const / (wavelength**2))

    # Total rotation in degrees
    total_rotation = specific_rotation * length_dm * concentration_g_ml

    print(f"For {color} light (Wavelength = {wavelength:.0f} nm):")
    # We were asked to output each number in the final equation.
    print(f"Total Rotation = ({A_const} + {B_const:.2e} / {wavelength**2:.0f}) * {length_dm} * {concentration_g_ml}")
    print(f"               = {specific_rotation:.2f} * {length_dm} * {concentration_g_ml}")
    print(f"               = {total_rotation:.1f} degrees")
    print("-" * 70)

print("The calculation shows that blue light is rotated significantly more than red light.")
print("This separation of colors along the tube's length, combined with the continuous")
print("rotation, results in a visual effect of a spiraling rainbow.")