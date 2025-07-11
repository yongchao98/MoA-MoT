# The problem asks for the cardinality of the set of colors used to construct the plane.
# Based on the logical deduction, the plane can only be constructed using White and Blue cubes.
# The set of colors is {White, Blue}.
# The cardinality of this set is the number of elements in it.

set_of_colors = {"White", "Blue"}
cardinality = len(set_of_colors)

# The final equation is simply the result of this count.
# We present the numbers that lead to the final answer.
# Let 1 represent White and 2 represent Blue.
color_1 = 1
color_2 = 2

# The total number of colors is the sum of the unique colors we have.
# In a set, this is just the count of items.
# Let's show it as an "equation" for clarity as requested.
print(f"The plane is built from 2 colors (White and Blue).")
print(f"Let's represent this as an equation of counting unique items:")
print(f"1 (for White) + 1 (for Blue) = {cardinality}")