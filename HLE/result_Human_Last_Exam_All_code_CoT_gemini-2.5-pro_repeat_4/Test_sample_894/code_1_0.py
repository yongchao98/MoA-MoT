# Propositions D and E are competing measures. If both pass, the one with fewer
# 'yes' votes is nullified. Thus, both could potentially be nullified.
# Propositions F and G are also competing measures with the same condition.
# Thus, both could also potentially be nullified.
propositions_that_could_be_nullified = ['D', 'E', 'F', 'G']

# The list is already in alphabetical order, but we sort it to be certain.
propositions_that_could_be_nullified.sort()

# Join the list into a comma-separated string with no spaces.
output_string = ",".join(propositions_that_could_be_nullified)

# Print the final result.
print(output_string)