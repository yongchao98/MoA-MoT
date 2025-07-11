import math

# Define the points based on the drawing instructions.
# I will use floating point numbers for consistency.
b1 = (0.0, -3.0)
b2 = (1.0, -3.0)
p = (2/3, -4.0)
c = (0.5, -1.5)
a1 = (2/3, -1.0)
a2 = (2/3, -3.0)

# Other points that define the shape
top_left = (0.0, 0.0)
top_right = (1.0, 0.0)
s_or_r = (1.0, -1.0) # The point on the top-right of the "rudder"

# The final drawing is made of several line segments.
# The "final equation" can be seen as the set of all numbers defining the points.
# The unique points are: b1, b2, p, c, a1, a2, top_left, top_right, s_or_r.

print("Numbers defining the final drawing's geometric equation:")

# Create a list of all unique points
points = {
    "top_left": top_left,
    "top_right": top_right,
    "s": s_or_r,
    "b1": b1,
    "p": p,
    "b2": b2,
    "a1": a1,
    "a2": a2,
    "c": c
}

# Print each number from the coordinates of the points
all_coords = []
for name, coords in points.items():
    # To avoid duplicates in the printout, we add numbers to a set first
    all_coords.append(coords[0])
    all_coords.append(coords[1])

# Using a set to keep only unique numbers, then sorting for clear output
unique_numbers = sorted(list(set(all_coords)))

# Output each unique number
for number in unique_numbers:
    print(f"{number:.4f}")
