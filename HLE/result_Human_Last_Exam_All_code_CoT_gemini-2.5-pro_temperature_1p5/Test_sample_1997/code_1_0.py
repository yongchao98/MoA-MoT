# The partition for Chapter 39 ("Yugiri") of the Tale of Genji is based on its incense pattern.
# The pattern connects sticks {1, 3} and {2, 4}, with {5} being separate.
# We represent this partition as a list of lists.
partition_data = [[1, 3], [2, 4], [5]]

# The data is already sorted as required (subsets sorted internally, and the list of subsets
# sorted by their first element). Now, we format it for printing.

# Build the string representation of the set of sets, e.g., "{{1,3},{2,4},{5}}"
inner_set_strings = []
for subset in partition_data:
    # Create strings for the inner sets, e.g., "{1,3}"
    subset_string = "{" + ",".join(map(str, subset)) + "}"
    inner_set_strings.append(subset_string)

# Join the inner set strings and wrap with the outer braces
final_output_string = "{" + ",".join(inner_set_strings) + "}"

print(final_output_string)