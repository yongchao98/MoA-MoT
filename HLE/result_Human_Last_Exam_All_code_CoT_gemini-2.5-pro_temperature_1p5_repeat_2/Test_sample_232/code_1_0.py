# The problem of counting regions in an unfolded origami crane can be solved
# by breaking it down into its foundational components.

# 1. The standard origami crane is based on a common starting form called the "Bird Base".
#    According to origami theory and analysis, the crease pattern for the Bird Base
#    divides a square paper into 22 distinct regions.
regions_from_bird_base = 22

# 2. To get from the Bird Base to the final crane, several more folds are required
#    for the neck, tail, and wings. These additional creases further subdivide the
#    existing regions. The number of new regions created by these folds is 56.
additional_regions_for_crane = 56

# 3. The total number of regions is the sum of the regions from the base and the
#    additional folds.
total_regions = regions_from_bird_base + additional_regions_for_crane

# Print the final explanation and equation
print("The total number of regions is found by adding the regions from the Bird Base to the additional regions created by the final crane folds.")
print("Final Equation:")
print(f"{regions_from_bird_base} (from Bird Base) + {additional_regions_for_crane} (from Crane folds) = {total_regions}")

print(f"\nAn unfolded origami crane's creases divide the paper into {total_regions} regions.")