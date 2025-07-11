# Number of fundamental segments connecting the two intersection points
# (line segment, upper arc, lower arc).
num_segments = 3

# The total number of subsets of these segments.
# This represents all possible combinations of segments a path could traverse.
total_subsets = 2**num_segments

# A valid path must traverse at least one segment, so we subtract 1
# to exclude the empty set (which corresponds to not traversing any segment).
num_distinct_paths = total_subsets - 1

# Print the final equation with each number explicitly shown.
print(f"The number of paths is based on the non-empty subsets of the {num_segments} available segments.")
print(f"The calculation is: 2 ** {num_segments} - 1 = {num_distinct_paths}")
print(f"So, there are {num_distinct_paths} distinct paths.")