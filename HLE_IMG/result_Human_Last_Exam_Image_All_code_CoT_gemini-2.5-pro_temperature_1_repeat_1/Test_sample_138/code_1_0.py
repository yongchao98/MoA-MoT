# The task is to identify the biological sexes of three pairs of insects (A, B, C)
# and provide the corresponding index from the given list of options.
#
# Options:
# 1) M, M
# 2) F, F
# 3) M, F
# 4) F, M
#
# Analysis:
# - Pair A (Cuckoo Bees): The left insect has a pointed abdomen (female ovipositor). The right insect has a pronged abdomen tip (male). So, A is F, M. This is option 4.
# - Pair B (Paper Wasps): The left insect has curled antennae (male). The right insect has straight antennae (female). So, B is M, F. This is option 3.
# - Pair C (Long-horned Bees): The left insect has very long antennae and a yellow face (male). The right insect has short antennae and a dark face (female). So, C is M, F. This is option 3.

# Assign the identified indices to variables.
index_A = 4
index_B = 3
index_C = 3

# Print the final result in the format "A, B, C".
print(f"{index_A}, {index_B}, {index_C}")