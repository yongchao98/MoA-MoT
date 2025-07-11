import numpy as np

# Step 1: Define the estimated infinity-norms for each of the 9 images.
# The norms are estimated by visually inspecting the images and finding the
# maximum absolute row sum. We assign relative magnitudes to the colors:
# Bright = 2.0, Medium = 1.0, Faint = 0.5.
# The sign of the value (red/green) is irrelevant for the infinity-norm.

# Estimated composition of the strongest row for each image:
# - Image 1: 4 bright, 4 medium => 4*2.0 + 4*1.0 = 12.0
# - Image 2: 4 bright, 2 medium, 2 faint => 4*2.0 + 2*1.0 + 2*0.5 = 11.0
# - Image 3: 4 bright, 4 medium => 4*2.0 + 4*1.0 = 12.0
# - Image 4: 4 bright, 2 medium, 2 faint => 4*2.0 + 2*1.0 + 2*0.5 = 11.0
# - Image 5: 5 bright, 2 medium, 2 faint => 5*2.0 + 2*1.0 + 2*0.5 = 13.0
# - Image 6: 2 bright, 2 medium, 2 faint => 2*2.0 + 2*1.0 + 2*0.5 = 7.0
# - Image 7: 3 bright, 2 medium, 2 faint => 3*2.0 + 2*1.0 + 2*0.5 = 9.0
# - Image 8: 3 bright, 2 medium, 2 faint => 3*2.0 + 2*1.0 + 2*0.5 = 9.0
# - Image 9: Appears as a lower-intensity version of Image 3. Assuming an
#            intensity scaling factor of 0.8, its norm is 12.0 * 0.8 = 9.6.

estimated_norms = [12.0, 11.0, 12.0, 11.0, 13.0, 7.0, 9.0, 9.0, 9.6]

# Step 2: Find the maximum and minimum norms from the list.
max_norm = max(estimated_norms)
min_norm = min(estimated_norms)

# Step 3: Calculate the ratio.
ratio = max_norm / min_norm

# Step 4: Print the final equation with each number and the result.
print("Based on visual estimation, the infinity-norms are:")
print(estimated_norms)
print(f"\nMaximum norm found: {max_norm}")
print(f"Minimum norm found: {min_norm}")
print("\nThe equation for the ratio is:")
print(f"Ratio = Maximum Norm / Minimum Norm")
print(f"Ratio = {max_norm} / {min_norm} = {ratio}")
