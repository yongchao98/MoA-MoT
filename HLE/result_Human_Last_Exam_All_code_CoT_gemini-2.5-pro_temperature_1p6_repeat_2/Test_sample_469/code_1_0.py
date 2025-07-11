import math

# Step 1: Define the properties of the system
ph = 4.21
concentration_mg_ml = 13
side_lengths_mm = [1, 100, 100] # converting 10cm to 100mm

# pH properties for Bromophenol blue
ph_range_low = 3.0
ph_range_high = 4.6
color_low_ph = "yellow"
color_high_ph = "blue"
color_mid_ph = "green"

# Step 2: Determine the intrinsic color of the solution based on pH
base_color = ""
if ph < ph_range_low:
    base_color = color_low_ph
elif ph > ph_range_high:
    base_color = color_high_ph
else:
    base_color = color_mid_ph

# Step 3: Define path lengths for thin and thick sides in cm
thin_path_cm = min(side_lengths_mm) / 10.0
thick_path_cm = max(side_lengths_mm) / 10.0

# Step 4: Explain the result using the Beer-Lambert Law principle
print(f"The given pH is {ph}.")
print(f"Bromophenol blue is {color_low_ph} below pH {ph_range_low} and {color_high_ph} above pH {ph_range_high}.")
print(f"At pH {ph}, it is in its transition range, so the solution's base color is {base_color}.")
print("-" * 20)
print("According to the Beer-Lambert Law, the intensity of observed color depends on the light's path length through the solution.")
print(f"Path length through the thin side = {thin_path_cm} cm.")
print(f"Path length through the thick side = {thick_path_cm} cm.")
print("-" * 20)
print(f"The path length through the thick side is {thick_path_cm / thin_path_cm:.0f} times longer than the thin side.")
print("A shorter path length results in a lighter, less intense color.")
print("A longer path length results in a darker, more intense color.")
print("-" * 20)
print(f"Conclusion for the thin side: The color will be light {base_color}.")
print(f"Conclusion for the thick side: The color will be a more intense {base_color}.")
print("\nThis corresponds to Answer Choice C.")
