# The Genji-kō pattern for Chapter 39 ("Yūgiri") connects the first, second,
# and fourth vertical lines (from the right), leaving the third and fifth
# lines independent.
# This corresponds to a partition of the set {1, 2, 3, 4, 5}.

# Define the subsets of the partition.
subset1 = {1, 2, 4}
subset2 = {3}
subset3 = {5}

# Create a list of the subsets.
partition = [subset1, subset2, subset3]

# Sort the list of subsets based on their minimum element to ensure
# they are in increasing order as requested.
# The lambda function x: min(x) gets the smallest number from each set.
partition.sort(key=lambda x: min(x))

# Prepare the output string.
# We iterate through the sorted partition and format each subset.
set_strings = []
for subset in partition:
    # Sort the numbers within each subset for consistent output, e.g., {1, 2, 4} not {2, 1, 4}
    subset_str = "{" + ",".join(map(str, sorted(list(subset)))) + "}"
    set_strings.append(subset_str)

# Join the formatted subset strings with commas and wrap them in outer braces.
final_output = "{" + ",".join(set_strings) + "}"

print(final_output)