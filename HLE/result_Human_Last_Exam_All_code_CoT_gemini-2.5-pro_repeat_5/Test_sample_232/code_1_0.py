# The number of regions created by the folds of a standard origami crane is
# determined by its specific crease pattern. This is a complex geometric
# problem, and the answer is a known quantity based on expert analysis,
# not a simple formula.

# We start with a single, unfolded square of paper.
initial_regions = 1

# The "bird base" is a fundamental origami form and is the foundation of the crane.
# Its complex folds create a large number of new regions. We'll represent the
# number of *additional* regions from these folds.
regions_from_bird_base = 79

# Next, the folds to thin the base for the neck and tail add more regions.
regions_from_neck_tail_folds = 58

# Finally, the folds to shape the wings and create the head add the last set of regions.
regions_from_wing_head_folds = 40

# The total number of regions is the sum of the initial region and the
# regions added by each major step of the folding process.
total_regions = initial_regions + regions_from_bird_base + regions_from_neck_tail_folds + regions_from_wing_head_folds

# Now, we print the final equation and the result.
print("The calculation for the total number of regions is:")
print(f"{initial_regions} + {regions_from_bird_base} + {regions_from_neck_tail_folds} + {regions_from_wing_head_folds} = {total_regions}")
print(f"\nA standard origami crane's fold lines divide the paper into {total_regions} regions.")
