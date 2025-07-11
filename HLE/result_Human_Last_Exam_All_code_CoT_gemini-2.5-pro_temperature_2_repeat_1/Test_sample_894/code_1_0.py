# This script identifies the San Francisco propositions from the November 2024 election
# that could be nullified by a competing proposition receiving more votes.

# A list of propositions with clauses that make them potentially nullifiable
# if their specific competing measure also passes but with a higher vote count.
# - C and D are competing measures.
# - E and F are competing measures.
# - J and K are competing measures.
propositions_that_could_be_nullified = ['C', 'D', 'E', 'F', 'J', 'K']

# Sort the list alphabetically to ensure correct order.
propositions_that_could_be_nullified.sort()

# Join the list of letters into a single string, separated by commas with no spaces.
output_string = ",".join(propositions_that_could_be_nullified)

# Print the final result.
print(output_string)