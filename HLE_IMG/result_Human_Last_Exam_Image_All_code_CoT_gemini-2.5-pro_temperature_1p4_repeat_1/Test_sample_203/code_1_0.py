# The Heesch number of a shape is the maximum number of layers of copies
# that can surround it. If a shape can tile the plane, we can add infinite
# layers, so its Heesch number is defined as infinity (∞).

# 1. Analyze the polygon:
# The image shows three identical polygons. We only need to analyze one.
# The polygon has point symmetry (180-degree rotational symmetry).
# - The V-shaped notch at the top is identical to the V-shaped notch at the bottom
#   but rotated 180 degrees.
# - The triangular protrusions on the top-left and bottom-right sides are a
#   180-degree rotated pair.
# - The triangular protrusions on the top-right and bottom-left sides are also
#   a 180-degree rotated pair.

# 2. Check for tiling ability:
# A known theorem in geometry states that any centrally symmetric polygon
# can tile the plane (form a tessellation).

# 3. Conclude the Heesch number:
# Since the polygon tiles the plane, its Heesch number is infinity (∞).

# 4. Format the final answer:
# The problem asks for the Heesch number for each of the three identical
# polygons, separated by commas.

heesch_number = "∞"
number_of_polygons = 3

# Create the final answer string
final_answer = ", ".join([heesch_number] * number_of_polygons)

print(final_answer)