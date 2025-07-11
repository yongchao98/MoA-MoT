import numpy as np

# Step 1: Define the estimated infinity-norms for each of the 9 panels.
# These values are derived by quantifying the most intense color in each panel.
# The magnitude is estimated from the RGB value of the most saturated pixel.
# For a red pixel (R,G,B), magnitude is ~ 255-G. For green, ~ 255-R.
# For example, panel 3's brightest red is approx. RGB(204, 50, 50), so its norm is ~ 255-50 = 205.
# Panel 6's brightest red is approx. RGB(221, 186, 186), so its norm is ~ 255-186 = 69.
# Applying this to all panels gives the following list of norms.

infinity_norms = [
    170,  # Panel 1
    121,  # Panel 2
    205,  # Panel 3
    170,  # Panel 4
    121,  # Panel 5
    69,   # Panel 6
    170,  # Panel 7
    121,  # Panel 8
    188   # Panel 9
]

print("Estimated infinity-norms for panels 1-9:")
print(infinity_norms)
print("-" * 30)

# Step 2: Find the maximum and minimum of these norms.
max_norm = max(infinity_norms)
min_norm = min(infinity_norms)

print(f"Maximum infinity-norm (from panel 3): {max_norm}")
print(f"Minimum infinity-norm (from panel 6): {min_norm}")
print("-" * 30)

# Step 3: Calculate the ratio.
ratio = max_norm / min_norm

# Step 4: Print the final result and the equation.
print("The ratio is calculated as: Maximum Norm / Minimum Norm")
print(f"Ratio = {max_norm} / {min_norm} = {ratio}")