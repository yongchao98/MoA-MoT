# Plan:
# 1. The set of visible trees that can act as blockers is {A, B, C, D}. The size of this set is 4.
# 2. A child's position is uniquely determined by an ordered pair of distinct blocking trees: one for tree E and one for tree F.
# 3. We calculate the number of ways to choose the tree that blocks the view to E. There are 4 choices.
# 4. We then calculate the number of ways to choose the tree that blocks the view to F. Since it must be different from the first, there are 3 remaining choices.
# 5. The maximum number of children is the product of these choices, which corresponds to the number of permutations of 2 items from a set of 4.

# Number of visible trees that can block the view
num_blockers = 4

# Number of choices for the tree blocking E
choices_for_blocking_E = num_blockers

# Number of choices for the tree blocking F (must be a different tree)
choices_for_blocking_F = num_blockers - 1

# The maximum number of children is the product of the choices
max_children = choices_for_blocking_E * choices_for_blocking_F

# The final equation is the number of choices for the first blocking tree
# multiplied by the number of choices for the second blocking tree.
# The code below prints each number in the final equation as requested.
print(f"{choices_for_blocking_E} * {choices_for_blocking_F} = {max_children}")