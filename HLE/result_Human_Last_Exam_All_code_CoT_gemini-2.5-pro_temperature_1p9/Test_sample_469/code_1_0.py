import math

# Step 1: Define the known values.
# The pH of the solution is given.
ph = 4.21
# The pKa of Bromophenol blue is a known constant.
pka_bromophenol_blue = 4.1
# The path lengths for the thin and thick sides.
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to cm for consistency.
path_length_thin_cm = path_length_thin_mm / 10.0

# The acid form (HIn) of Bromophenol blue is yellow.
# The base form (In-) of Bromophenol blue is blue.

# Step 2: Use the Henderson-Hasselbalch equation to find the ratio of the base (blue) form to the acid (yellow) form.
# The equation is: pH = pKa + log10([Base]/[Acid])
# Rearranging for the ratio: [Base]/[Acid] = 10^(pH - pKa)
ratio_blue_to_yellow = 10**(ph - pka_bromophenol_blue)

# Step 3: Print the results and explain the implications.
print("--- Step 1: Determine Solution Color ---")
print(f"The equation for the ratio of indicator forms is: [Blue Form] / [Yellow Form] = 10^(pH - pKa)")
# Remember to output each number in the final equation!
print(f"Plugging in the values: [Blue Form] / [Yellow Form] = 10^({ph} - {pka_bromophenol_blue})")
print(f"Calculated ratio of [Blue Form] / [Yellow Form]: {ratio_blue_to_yellow:.2f}")
print("\nSince the ratio is close to 1, both the yellow and blue forms of the indicator are present in significant amounts.")
print("The mixture of blue and yellow produces a GREEN solution.")

print("\n--- Step 2: Determine Color Intensity ---")
print("According to the Beer-Lambert Law, color intensity is proportional to the path length of light.")
print(f"Path length of the thin side: {path_length_thin_cm} cm")
print(f"Path length of the thick side: {path_length_thick_cm} cm")
print("The path length of the thick side is 100 times greater than the thin side.")
print("This will result in a much more intense color when viewed through the thick side.")

print("\n--- Conclusion ---")
print("Thin side: The short path length results in a less intense color, which will appear as LIGHT GREEN.")
print("Thick side: The long path length results in a more intense color, which will appear as GREEN.")
