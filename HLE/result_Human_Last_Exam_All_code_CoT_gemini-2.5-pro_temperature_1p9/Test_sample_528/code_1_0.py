# This script calculates the number of sets (S1, S2, S3, S4) satisfying the given conditions.
# The method is based on counting the possible placements for each element
# in the universal set U = {1, 2, 3, 4, 5}.

# The chain of subsets S1
# creates 5 disjoint regions. We determine the number of choices for placing each element
# of U into these regions.

# Number of choices for placing element 1, given the condition 1 in S1.
# This means 1 must be in S1.
choices_for_1 = 1

# Number of choices for placing element 2, given the condition 2 in S2.
# This means 2 can be in S1 or (S2 \ S1).
choices_for_2 = 2

# Number of choices for placing element 3, given the condition 3 in S3.
# This means 3 can be in S1, (S2 \ S1), or (S3 \ S2).
choices_for_3 = 3

# Number of choices for placing element 4. There are no restrictions.
# It can be in any of the 5 regions: S1, S2\S1, S3\S2, S4\S3, U\S4.
choices_for_4 = 5

# Number of choices for placing element 5. There are no restrictions.
# It can also be in any of the 5 regions.
choices_for_5 = 5

# The total number of configurations is the product of the choices for each element.
total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

print("The number of choices for placing each element is as follows:")
print("Element 1: {}".format(choices_for_1))
print("Element 2: {}".format(choices_for_2))
print("Element 3: {}".format(choices_for_3))
print("Element 4: {}".format(choices_for_4))
print("Element 5: {}".format(choices_for_5))
print("\nThe total number of sets is the product of these choices:")
print("Total number = {} * {} * {} * {} * {}".format(choices_for_1, choices_for_2, choices_for_3, choices_for_4, choices_for_5))
print("Result = {}".format(total_sets))