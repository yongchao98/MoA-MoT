import math

# Step 1: Define the known values.
pKa_bromophenol_blue = 4.1
pH_solution = 4.21
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to cm for consistency
path_length_thin_cm = path_length_thin_mm / 10.0
path_length_thick_cm = 10.0

print("Step-by-step Analysis:")
print("-----------------------")
print(f"The pH of the solution is {pH_solution}.")
print(f"The pKa of Bromophenol blue is approximately {pKa_bromophenol_blue}.")
print("\nThe color of a pH indicator depends on the ratio of its acidic and basic forms.")
print("The acidic form (HIn) of Bromophenol blue is yellow.")
print("The basic form (In-) of Bromophenol blue is blue.")

# Step 2: Use the Henderson-Hasselbalch equation to find the ratio of the forms.
# pH = pKa + log10([In-]/[HIn])
# log10([In-]/[HIn]) = pH - pKa
# [In-]/[HIn] = 10^(pH - pKa)

ratio_log = pH_solution - pKa_bromophenol_blue
ratio = 10**ratio_log

print("\nUsing the Henderson-Hasselbalch equation:")
print(f"log([Blue Form]/[Yellow Form]) = pH - pKa = {pH_solution} - {pKa_bromophenol_blue} = {ratio_log:.2f}")
print(f"The ratio of [Blue Form]/[Yellow Form] = 10^{ratio_log:.2f} = {ratio:.2f}")
print("\nBecause the ratio is close to 1, there are significant amounts of both the yellow and blue forms.")
print("The mixture of yellow and blue light appears green to the human eye.")

# Step 3: Apply the Beer-Lambert Law to consider the path length.
print("\nAccording to the Beer-Lambert Law (Absorbance = constant * path_length * concentration), the intensity of the color depends on the path length of the light.")

print(f"\nCase 1: Viewing through the thin side.")
print(f"Path length = {path_length_thin_mm} mm ({path_length_thin_cm} cm).")
print("With a short path length, the absorbance is low. The color will be a light green.")

print(f"\nCase 2: Viewing through the thick side.")
print(f"Path length = {path_length_thick_cm} cm.")
print(f"This path length is {path_length_thick_cm / path_length_thin_cm} times longer than the thin side.")
print("With a long path length, the absorbance is much higher. The color will be a much more intense, darker green.")

print("\n-----------------------")
print("Conclusion:")
print("Thin side color: Light Green")
print("Thick side color: Green")
<<<C>>>