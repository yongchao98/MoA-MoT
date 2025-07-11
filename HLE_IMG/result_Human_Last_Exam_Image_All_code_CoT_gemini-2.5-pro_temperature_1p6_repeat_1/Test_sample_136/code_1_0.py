# The number of beads on the necklace defines the order of the group.
n_beads = 6

# The group of symmetries is identified as the cyclic group C6 because each equivalence
# class (row) contains only rotated versions of a single coloring. Reflections are excluded.
# For a cyclic group of order n, a minimal generator is the rotation by the smallest angle.
# The calculation for this angle is: 360 / n
rotation_angle = 360 / n_beads

# The problem asks for a comma-separated list of the minimal generators.
# For the cyclic group C6, there is only one minimal generator required.
# (Note: rotation by 300 degrees, or -60 degrees, is also a valid generator,
# but by convention, we choose the smallest positive angle).
generator_description = f"rotation by {int(rotation_angle)} degrees"

print("The minimal generating set for the group of symmetries is:")
print(generator_description)