# The incense pattern for Chapter 39, "Yūgiri", of the Tale of Genji
# corresponds to a specific partition of the set {1, 2, 3, 4, 5}.
#
# The pattern connects the first two sticks and the third and fourth sticks,
# leaving the fifth stick independent.
#
# This creates the following groups:
# - Group 1: {1, 2}
# - Group 2: {3, 4}
# - Group 3: {5}
#
# We will represent this as a list of lists, sorted as required.
partition_chap39 = [[1, 2], [3, 4], [5]]

# Format the partition into the required string representation "{{1, 2}, {3, 4}, {5}}".
inner_sets_str = []
for subset in partition_chap39:
    # Format each inner list like "{n1, n2, ...}"
    inner_str = "{" + ", ".join(map(str, subset)) + "}"
    inner_sets_str.append(inner_str)

# Combine the inner set strings and wrap them in the outer braces.
final_output_str = "{" + ", ".join(inner_sets_str) + "}"

print(final_output_str)