import math

# Define the known parameters
pKa_bromophenol_blue = 4.1
pH_solution = 4.21
path_length_thin_mm = 1
path_length_thick_cm = 10

# Convert all path lengths to cm for consistency
path_length_thin_cm = path_length_thin_mm / 10.0

print("This script determines the color of a Bromophenol blue solution based on pH and viewing depth.")
print("-" * 70)

# --- Step 1: Determine the fundamental color of the solution ---
print("Step 1: Determine the color of the Bromophenol Blue solution at the given pH.")
print("The color of a pH indicator depends on the ratio of its acidic and basic forms.")
print("This ratio can be calculated using the Henderson-Hasselbalch equation: pH = pKa + log([Base]/[Acid])")
print("\nThe acidic form of Bromophenol Blue is yellow.")
print("The basic form of Bromophenol Blue is blue.")
print(f"Given pH = {pH_solution} and the indicator's pKa = {pKa_bromophenol_blue}.")

# Calculate the ratio [Base]/[Acid]
ratio = 10**(pH_solution - pKa_bromophenol_blue)

print(f"\nThe ratio of [Base]/[Acid] is 10^({pH_solution} - {pKa_bromophenol_blue}) = {ratio:.2f}")
print("Since the ratio is close to 1, both the yellow (acidic) and blue (basic) forms are present.")
print("A mixture of yellow and blue light results in a perceived color of GREEN.")
print("-" * 70)

# --- Step 2: Analyze the effect of path length ---
print("Step 2: Compare the color intensity when viewed through the thin and thick sides.")
print("The Beer-Lambert Law states that absorbance is directly proportional to the path length of light.")
print(f"Path length of the thin side = {path_length_thin_cm} cm")
print(f"Path length of the thick side = {path_length_thick_cm} cm")
print("\nBecause the path length is much longer through the thick side, the absorbance will be much higher.")
print("Higher absorbance means the color appears more intense and saturated (less light passes through).")
print("-" * 70)

# --- Step 3: Conclusion ---
print("Conclusion:")
print("The fundamental color of the solution is GREEN.")
print("- Through the thin side (0.1 cm), the color will be a less intense green, described as 'light green'.")
print("- Through the thick side (10 cm), the color will be a much more intense and saturated green, described as 'green' or 'dark green'.")
print("\nBased on this analysis, the best answer choice is C.")

print("\n<<<C>>>")