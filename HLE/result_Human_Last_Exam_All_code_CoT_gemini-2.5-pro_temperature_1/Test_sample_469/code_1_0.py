import math

# Step 1: Define the known parameters for the chemical properties.
pKa_bromophenol_blue = 4.1  # pKa of Bromophenol blue
pH_solution = 4.21

# Step 2: Define the physical dimensions of the prism.
path_length_thin_mm = 1  # in mm
path_length_thick_cm = 10 # in cm

# Convert all path lengths to the same unit (cm) for comparison.
path_length_thin_cm = path_length_thin_mm / 10.0

# Step 3: Determine the intrinsic color of the solution using the Henderson-Hasselbalch equation.
# The equation relates pH, pKa, and the ratio of the base (blue) and acid (yellow) forms.
# Equation: pH = pKa + log10([Base]/[Acid])
# Rearranging to solve for the ratio: [Base]/[Acid] = 10^(pH - pKa)
ratio_base_to_acid = 10**(pH_solution - pKa_bromophenol_blue)

print("--- Part 1: Determining the Intrinsic Color ---")
print(f"The pH of the solution is {pH_solution} and the pKa of the Bromophenol blue indicator is {pKa_bromophenol_blue}.")
print("Using the Henderson-Hasselbalch equation, we calculate the ratio of the blue (base) form to the yellow (acid) form.")
print(f"Ratio = 10^({pH_solution} - {pKa_bromophenol_blue}) = {ratio_base_to_acid:.2f}")
print("Since the ratio of the blue form to the yellow form is close to 1, both colors are present in significant amounts.")
print("The combination of yellow and blue light makes the solution appear green.\n")


# Step 4: Compare the color intensity based on path length using the Beer-Lambert Law.
# The Beer-Lambert Law states that Absorbance is directly proportional to path length (l).
# A higher absorbance results in a more intense color.
absorbance_ratio = path_length_thick_cm / path_length_thin_cm

print("--- Part 2: Comparing Color Intensity by Path Length ---")
print(f"The path length of the thin side is {path_length_thin_mm} mm, which is {path_length_thin_cm} cm.")
print(f"The path length of the thick side is {path_length_thick_cm} cm.")
print("According to the Beer-Lambert Law, the intensity of the observed color is proportional to the path length.")
print(f"The path length through the thick side is {path_length_thick_cm} / {path_length_thin_cm} = {absorbance_ratio:.0f} times longer than the thin side.")
print("This means the color will appear much more intense (darker) when viewed through the thick side.\n")

# Step 5: Conclusion
print("--- Conclusion ---")
print("The solution has an intrinsic green color.")
print(f"Viewed through the thin side ({path_length_thin_cm} cm), the color will be light green due to low absorbance.")
print(f"Viewed through the thick side ({path_length_thick_cm} cm), the color will be a more intense green due to high absorbance.")
print("This matches the choice: Thin: light green, Thick: green.")

<<<C>>>