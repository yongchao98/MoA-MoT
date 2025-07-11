import math

# --- Setup ---
# Define the known values for the problem.
pKa_bromophenol_blue = 4.0
pH_solution = 4.21
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to a consistent unit (cm).
path_length_thin_cm = path_length_thin_mm / 10.0

# --- Analysis ---
print("Step 1: Determine the intrinsic color of the solution using the Henderson-Hasselbalch equation.")
print("The acidic form of Bromophenol blue is yellow, and the basic form is blue.")
print(f"The pKa is ~{pKa_bromophenol_blue} and the solution pH is {pH_solution}.")

# Calculate the ratio of the basic (blue) form to the acidic (yellow) form.
# Equation: pH = pKa + log10([Base]/[Acid])
log_ratio = pH_solution - pKa_bromophenol_blue
ratio_base_to_acid = 10**log_ratio

print("\nUsing the Henderson-Hasselbalch Equation:")
print(f"  {pH_solution} = {pKa_bromophenol_blue} + log([Blue Form]/[Yellow Form])")
print(f"  log([Blue Form]/[Yellow Form]) = {pH_solution} - {pKa_bromophenol_blue} = {log_ratio:.2f}")
print(f"  Ratio of [Blue Form]/[Yellow Form] = 10^{log_ratio:.2f} = {ratio_base_to_acid:.2f}")

print("\nConclusion for Step 1:")
print("Since both yellow and blue forms are present, their colors mix to appear green.")

print("\nStep 2: Determine the effect of viewing path length using the Beer-Lambert Law.")
print("The Beer-Lambert Law states that absorbance is proportional to path length.")
print(f"Path length of the thin side: {path_length_thin_cm} cm")
print(f"Path length of the thick side: {path_length_thick_cm} cm")

print("\nConclusion for Step 2:")
print("A shorter path length (the thin side) means lower absorbance, resulting in a light green color.")
print("A longer path length (the thick side) means higher absorbance, resulting in a more intense, darker green color.")

print("\n--- Final Answer ---")
print("Thin side color: Light Green")
print("Thick side color: Green")