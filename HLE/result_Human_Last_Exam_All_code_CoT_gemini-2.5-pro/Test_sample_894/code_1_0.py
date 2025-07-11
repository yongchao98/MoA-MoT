# This script identifies and lists the San Francisco propositions from the
# November 2024 election that would be nullified if a competing measure
# passes with more votes.

# Based on the official voter guide, the following pairs of propositions
# are in direct conflict:
# 1. Prop C vs. Prop D (Real Estate Transfer Tax)
# 2. Prop E vs. Prop G (Police Policies)
# 3. Prop I vs. Prop J (JFK Drive Access)

# In each case, if both propositions pass, the one with fewer votes is nullified.
# Therefore, any proposition in these pairs could potentially be nullified.
propositions_that_could_be_nullified = ['C', 'D', 'E', 'G', 'I', 'J']

# Sort the list alphabetically to ensure correct order.
propositions_that_could_be_nullified.sort()

# Join the list into a comma-separated string with no spaces.
final_list_string = ",".join(propositions_that_could_be_nullified)

# Print the final formatted list.
print(final_list_string)