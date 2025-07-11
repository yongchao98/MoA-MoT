# The rotational symmetry of a pattern is determined by finding a point
# around which the pattern can be rotated by a certain angle (less than 360 degrees)
# to look identical to its original state. The number of such positions in a full
# 360-degree rotation is the order of the rotational symmetry.

# By observing the provided image, we can identify central motifs that look like stars or rosettes.
# Let's count the number of repeating units (or "arms") in one of these rosettes.
# The rosette is composed of a cluster of orange rhombi all pointing to a central point.

# Counting the number of these orange rhombi gives us the order of rotation.
number_of_repeating_units = 7

# This means the pattern is unchanged after a rotation of 360 / 7 degrees.
# The order of rotational symmetry is therefore 7.

print(f"The central motif in the tiling has {number_of_repeating_units} repeating units.")
print(f"Therefore, the rotational symmetry of the tiling is {number_of_repeating_units}.")
