# This script calculates the total number of equivalence classes of endomorphisms on a set of size 4.
# The calculation is based on partitioning the 4 elements into components of a functional graph.

# Number of structures for a single component of size 4
count_part_4 = 9

# Number of structures for components of sizes 3 and 1
count_part_3_1 = 4

# Number of structures for two components of size 2
count_part_2_2 = 3

# Number of structures for components of sizes 2, 1, and 1
count_part_2_1_1 = 2

# Number of structures for four components of size 1
count_part_1_1_1_1 = 1

# Summing the counts from all partitions
total_classes = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1

# The problem asks to output the equation with each number.
print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_classes}")