# The plan is to determine the index for each pair based on morphological differences.
# Pair A: The left insect has a pointed abdomen (female), the right one has a blunt tip (male). This is F, M -> index 4.
# Pair B: The left insect has straight antennae (female), the right one has curled antennae (male). This is F, M -> index 4.
# Pair C: The left insect has very long antennae (male), the right one has short antennae (female). This is M, F -> index 3.

# Assigning the identified indices to variables.
answer_A = 4
answer_B = 4
answer_C = 3

# Printing the final answer in the requested format.
print(f"{answer_A}, {answer_B}, {answer_C}")