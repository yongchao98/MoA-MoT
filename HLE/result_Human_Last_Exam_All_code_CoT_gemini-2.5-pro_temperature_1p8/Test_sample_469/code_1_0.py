import math

# Step 1: Define the properties and given values.
# The pKa is the pH at which the acidic and basic forms of the indicator are in equal concentration.
pka_bromophenol_blue = 4.1
solution_ph = 4.21

# Path lengths of the rectangular prism.
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to cm for consistency.
path_length_thin_cm = path_length_thin_mm / 10.0

print(f"The pH of the solution is {solution_ph}.")
print(f"The pKa of Bromophenol Blue is {pka_bromophenol_blue}.")
print("-" * 20)
print("Step 1: Determine the intrinsic color of the solution.")

# Step 2: Use the Henderson-Hasselbalch equation to find the ratio of the basic (blue) to acidic (yellow) forms.
# pH = pKa + log10([Blue Form] / [Yellow Form])
# so, [Blue Form] / [Yellow Form] = 10^(pH - pKa)
ratio_blue_to_yellow = 10**(solution_ph - pka_bromophenol_blue)

print("The Henderson-Hasselbalch equation is: pH = pKa + log10([Blue Form]/[Yellow Form])")
print(f"Calculating the ratio: 10^({solution_ph} - {pka_bromophenol_blue}) = {ratio_blue_to_yellow:.2f}")
print(f"The ratio of the blue form to the yellow form is approximately {ratio_blue_to_yellow:.2f} to 1.")
print("Since both the yellow and blue forms are present, the solution's intrinsic color is green.")
print("-" * 20)

# Step 3: Analyze the effect of path length using the Beer-Lambert Law (A = ebc).
# Absorbance (A) is directly proportional to path length (b).
print("Step 2: Analyze the effect of path length on color intensity.")
print("According to the Beer-Lambert law, absorbance is proportional to the path length of light.")
print(f"Path length of the thin side: {path_length_thin_cm} cm")
print(f"Path length of the thick side: {path_length_thick_cm} cm")

# Calculate the ratio of the path lengths
path_ratio = path_length_thick_cm / path_length_thin_cm
print(f"\nThe absorbance through the thick side is {path_ratio:.0f} times greater than through the thin side.")
print("Higher absorbance results in a darker, more saturated color.")
print("-" * 20)

# Step 4: Conclude the final perceived colors.
print("Conclusion:")
print(f"Through the thin side ({path_length_thin_cm} cm), the absorbance is lower, so the color will be a light green.")
print(f"Through the thick side ({path_length_thick_cm} cm), the absorbance is much higher, so the color will be a saturated, dark green.")
print("\nComparing these results to the answer choices, 'Thin: light green, Thick: green' is the best description.")

<<<C>>>