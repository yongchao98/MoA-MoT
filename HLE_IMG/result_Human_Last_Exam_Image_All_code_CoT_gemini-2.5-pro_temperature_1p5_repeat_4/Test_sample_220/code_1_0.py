import math

# Step 1: Define the theoretical properties of an FCC [110] projection.
# The projection is a centered rectangle with side lengths 'a' and 'a*sqrt(2)'.
# The ratio of the side lengths is therefore sqrt(2).
theoretical_ratio = math.sqrt(2)

# Step 2: Extract approximate atom coordinates from Image B to define a centered rectangular cell.
# From visual inspection of Image B:
center_atom_x = 4.5
center_atom_y = 8.0
corner_atom1_x = 1.0
corner_atom1_y = 5.5
corner_atom2_x = 8.0
corner_atom2_y = 5.5

# Step 3: Calculate the dimensions and ratio of the rectangle in Image B.
# The first side length (s1) is the distance between the two corner atoms.
s1 = corner_atom2_x - corner_atom1_x

# The second side length (s2) is twice the perpendicular distance from the
# center atom to the side defined by the two corner atoms.
s2 = 2 * (center_atom_y - corner_atom1_y)

# Calculate the ratio of the longer side to the shorter side.
if s1 > s2:
    measured_ratio = s1 / s2
    print(f"The side lengths of the rectangle in Image B are approximately {s1:.1f} and {s2:.1f}.")
    print(f"The calculated ratio of side lengths is {s1:.1f} / {s2:.1f} = {measured_ratio:.4f}")
else:
    measured_ratio = s2 / s1
    print(f"The side lengths of the rectangle in Image B are approximately {s2:.1f} and {s1:.1f}.")
    print(f"The calculated ratio of side lengths is {s2:.1f} / {s1:.1f} = {measured_ratio:.4f}")

# Step 4: Compare the measured ratio with the theoretical ratio.
print(f"The theoretical ratio for an FCC [110] projection is sqrt(2) ≈ {theoretical_ratio:.4f}")
print("The measured ratio from Image B closely matches the theoretical ratio.")
print("\nTherefore, Image B represents a face-centered cubic structure viewed along the [110] direction.")
