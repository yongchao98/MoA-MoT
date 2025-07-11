import math

# Step 1: Define the problem's constants.
# The pH of the aqueous Bromophenol blue solution.
ph_solution = 4.21
# The pKa of Bromophenol blue, around which it changes color, is approximately 4.1.
pka_indicator = 4.1
# The dimensions of the rectangular prism give the path lengths.
# "Thin side" path length is 1 mm, which is 0.1 cm.
path_length_thin = 0.1
# "Thick side" path length is 10 cm.
path_length_thick = 10.0

# Step 2: Use the Henderson-Hasselbalch equation to determine the intrinsic color.
# The equation is: pH = pKa + log10([Basic Form] / [Acidic Form])
# Bromophenol blue is yellow in its acidic form and blue in its basic form.
# We can calculate the ratio of the blue form to the yellow form.
ratio_blue_to_yellow = 10**(ph_solution - pka_indicator)

print("--- Analysis of Solution Color ---")
print("The Henderson-Hasselbalch equation is used to find the ratio of the basic (blue) form to the acidic (yellow) form of the indicator.")
print(f"Ratio = 10^(pH - pKa)")
print(f"Ratio = 10^({ph_solution} - {pka_indicator}) = {ratio_blue_to_yellow:.2f}")
print("\nBecause the ratio is close to 1, the solution contains a significant mixture of both the yellow and blue forms.")
print("A mixture of yellow and blue light appears GREEN. So, the intrinsic color of the solution is green.")

# Step 3: Apply the Beer-Lambert Law to explain the effect of path length.
# Beer-Lambert Law: Absorbance = ε * b * c (where b is path length)
# A longer path length results in higher absorbance, making the color appear more intense/darker.
# A shorter path length results in lower absorbance, making the color appear lighter.
print("\n--- Analysis of Perceived Color vs. Path Length ---")
print(f"Path length of the 'thin side': {path_length_thin} cm")
print(f"Path length of the 'thick side': {path_length_thick} cm")
print(f"The thick side's path is {path_length_thick / path_length_thin} times longer than the thin side's.")
print("\nAccording to the Beer-Lambert law, a longer path length leads to higher light absorbance.")
print("- Through the thin side (short path), the absorbance is low. The color will appear LIGHT GREEN.")
print("- Through the thick side (long path), the absorbance is high. The color will appear as a more saturated and darker GREEN.")

print("\n--- Final Answer ---")
print("Thin: light green, Thick: green")