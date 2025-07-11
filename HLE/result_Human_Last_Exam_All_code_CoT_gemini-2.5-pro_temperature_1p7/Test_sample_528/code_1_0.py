# The problem is to find the number of set tuples (S1, S2, S3, S4) such that
# S1 ⊂ S2 ⊂ S3 ⊂ S4 ⊂ {1, 2, 3, 4, 5} and i ∈ Si for i = 1, 2, 3.

# We determine the number of choices for placing each element from {1, 2, 3, 4, 5}
# into the chain of sets.

# For element 1, the condition is 1 ∈ S1. There is only 1 choice.
choices_for_1 = 1

# For element 2, the condition is 2 ∈ S2. Element 2 can be in S1 or (S2 \ S1). There are 2 choices.
choices_for_2 = 2

# For element 3, the condition is 3 ∈ S3. Element 3 can be in S1, (S2 \ S1), or (S3 \ S2). There are 3 choices.
choices_for_3 = 3

# For element 4, there are no constraints. It can be in S1, (S2 \ S1), (S3 \ S2), (S4 \ S3), or outside S4. There are 5 choices.
choices_for_4 = 5

# For element 5, there are no constraints. It also has 5 choices.
choices_for_5 = 5

# The total number of valid set tuples is the product of the choices for each element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

# Print the final equation as requested.
print(f"The number of sets is the product of the number of choices for each element.")
print(f"Number of choices for each element are: {choices_for_1}, {choices_for_2}, {choices_for_3}, {choices_for_4}, {choices_for_5}")
print(f"The final calculation is:")
print(f"{choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")