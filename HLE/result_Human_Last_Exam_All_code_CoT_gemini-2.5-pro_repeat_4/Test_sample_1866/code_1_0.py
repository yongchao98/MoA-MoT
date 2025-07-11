import math

# Problem parameters
cube_size = 4  # The side length of the initial cube in cm
target_size = 1  # The side length of the final cubes in cm
knife_depth = 2  # The maximum cutting depth of the knife in cm

# --- Stage 1: Cut the 4x4x4 cube into 2x2x2 cubes ---
print("Stage 1: Cutting the 4x4x4 cube into 2x2x2 cubes.")

# To do this, we must make a cut through the center of each of the 3 dimensions.
# The thickness of the material for these first cuts is the initial cube size.
thickness_stage1 = cube_size

# Since the knife can only cut 2cm deep, we calculate the number of cuts needed
# to get through the 4cm thickness. We can flip the cube, so we need 2 cuts.
cuts_per_central_plane = math.ceil(thickness_stage1 / knife_depth)
num_dimensions = 3
cuts_stage1 = num_dimensions * cuts_per_central_plane

print(f"The initial cube is {thickness_stage1}x{thickness_stage1}x{thickness_stage1} cm.")
print(f"The knife can only cut {knife_depth} cm deep.")
print(f"To cut through a {thickness_stage1} cm block, we need {cuts_per_central_plane} cuts (one from each side).")
print(f"This must be done for all {num_dimensions} dimensions (X, Y, Z).")
print(f"Cuts for Stage 1 = {num_dimensions} * {cuts_per_central_plane} = {cuts_stage1}")
print("-" * 20)

# After stage 1, we have 8 cubes of size 2x2x2.

# --- Stage 2: Cut the 2x2x2 cubes into 1x1x1 cubes ---
print("Stage 2: Cutting the 2x2x2 cubes into 1x1x1 cubes.")

# Now we need to cut each 2x2x2 cube in half along each of its 3 axes.
# The thickness of the material for these cuts is 2cm.
thickness_stage2 = cube_size / cuts_per_central_plane # This is 2cm

# The knife can cut 2cm in a single pass.
# By arranging all the pieces side-by-side, we can make all the cuts
# for one dimension in a single operation.
cuts_per_remaining_dimension = 1
cuts_stage2 = num_dimensions * cuts_per_remaining_dimension

print(f"After Stage 1, we have 8 cubes of size {int(thickness_stage2)}x{int(thickness_stage2)}x{int(thickness_stage2)} cm.")
print("To get 1x1x1 cubes, we must cut each of these in half along each of the 3 axes.")
print("Since the thickness is 2cm, the knife can cut through in a single pass.")
print("We can arrange all pieces to perform all cuts for one dimension in a single go.")
print(f"Cuts for Stage 2 = {num_dimensions} dimensions * {cuts_per_remaining_dimension} cut each = {cuts_stage2}")
print("-" * 20)

# --- Final Calculation ---
total_cuts = cuts_stage1 + cuts_stage2
print("The minimum total number of cuts is the sum of cuts from both stages.")
print(f"Total cuts = {cuts_stage1} (Stage 1) + {cuts_stage2} (Stage 2) = {total_cuts}")
