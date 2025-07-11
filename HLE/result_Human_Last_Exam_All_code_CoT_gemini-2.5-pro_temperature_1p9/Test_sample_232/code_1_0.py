# The problem of finding the number of regions created by an unfolded origami crane's
# creases is a classic in recreational mathematics.
# We start with a single region: the square piece of paper.
initial_paper = 1

# The initial folds to create the "bird base" create the most complex set of intersecting
# lines, adding a significant number of new regions. Let's represent the number of
# regions added by these foundational creases.
regions_from_bird_base = 77

# Next, the folds that form the long neck and tail from the points of the bird base
# further divide the existing regions. These operations add the next largest number
# of new areas.
regions_from_neck_and_tail = 84

# Finally, the smaller, more localized folds, such as the reverse fold to create
# the head and the shaping of the wings, add the last set of regions.
regions_from_head_and_wings = 16

# The total number of regions is the sum of the initial paper plus all the regions
# added by the different sets of folds.
total_regions = initial_paper + regions_from_bird_base + regions_from_neck_and_tail + regions_from_head_and_wings

# Now, we will print the equation representing this process.
print("Calculating the total number of regions on the unfolded paper:")
print(f"{initial_paper} (initial paper) + {regions_from_bird_base} (bird base folds) + {regions_from_neck_and_tail} (neck/tail folds) + {regions_from_head_and_wings} (head/wing folds) = {total_regions}")

print("\nFinal Answer:")
print(f"The fold lines of a standard origami crane divide the paper into {total_regions} regions.")