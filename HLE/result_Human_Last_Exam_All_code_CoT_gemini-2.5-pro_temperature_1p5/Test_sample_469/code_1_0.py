import math

# Step 1: Define the properties of the system
pKa_bromophenol_blue = 4.1
pH_solution = 4.21
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to the same unit (cm) for comparison
path_length_thin_cm = path_length_thin_mm / 10.0

# Step 2: Determine the ratio of indicator forms using the Henderson-Hasselbalch equation
# pH = pKa + log10([Basic Form] / [Acidic Form])
# The acidic form of Bromophenol blue is yellow, and the basic form is blue.
ratio_blue_to_yellow = 10**(pH_solution - pKa_bromophenol_blue)

# Step 3: Analyze the color based on the ratio of forms
# A ratio close to 1 means significant amounts of both forms are present.
# The mixing of yellow and blue light is perceived as green.
print("Step 1: Determine the base color of the solution.")
print(f"The pH of the solution is {pH_solution}.")
print(f"The pKa of Bromophenol blue is {pKa_bromophenol_blue}.")
print(f"The ratio of the blue form to the yellow form is 10^({pH_solution} - {pKa_bromophenol_blue}) = {ratio_blue_to_yellow:.2f}.")
print("Because both yellow and blue forms are present, the solution's fundamental color is green.\n")

# Step 4: Analyze the effect of path length on color intensity using the Beer-Lambert Law (A = ebc)
# Absorbance (A) is directly proportional to path length (b).
print("Step 2: Compare the color intensity based on path length.")
print(f"Path length through the thin side: {path_length_thin_cm} cm")
print(f"Path length through the thick side: {path_length_thick_cm} cm")

# Calculate the ratio of absorbance
absorbance_ratio = path_length_thick_cm / path_length_thin_cm
print(f"The path length through the thick side is {absorbance_ratio:.0f} times longer than through the thin side.")
print("According to the Beer-Lambert law, a longer path length results in higher absorbance and a more intense color.\n")

# Step 5: Conclude the final perceived colors
print("Step 3: Conclusion")
print("The base color of the solution is green.")
print("Through the thin side (short path length), the color will be faint: light green.")
print("Through the thick side (long path length), the color will be intense: green.")
print("\nTherefore, the correct answer is Thin: light green, Thick: green.")
