# This script determines the partition for Chapter 39 of The Tale of Genji.
# The incense pattern for Chapter 39 ("Yugiri") has the first stick separate
# and the remaining four sticks (2, 3, 4, 5) grouped together.

# We represent this partition as a list of lists.
partition_data = [[1], [2, 3, 4, 5]]

# For consistency and adherence to the specified format, we ensure all subsets
# and the main partition are sorted.
for subset in partition_data:
    subset.sort()
partition_data.sort()

# Now, we format this data into the required string representation {{...},{...}}.
# This loop processes each number from the partition data to build the final string.
subset_strings = []
for subset in partition_data:
    # Convert each number in the subset to a string and join with commas
    subset_content = ",".join(map(str, subset))
    # Enclose the content in braces to form a set string like {1} or {2,3,4,5}
    subset_strings.append("{" + subset_content + "}")

# Join the formatted subset strings with commas
final_string = "{" + ",".join(subset_strings) + "}"

print(final_string)