# This script identifies the San Francisco propositions on the November 2024 ballot
# that could be nullified by a competing measure.

# The pairs of competing measures are:
# - Prop C vs. Prop D
# - Prop J vs. Prop L
# - Prop M vs. Prop N
# In each case, if both pass, the one with fewer 'yes' votes is nullified.
# Therefore, all six are on the list of propositions that could be nullified.

nullifiable_propositions = ['C', 'D', 'J', 'L', 'M', 'N']

# Sort the list alphabetically
nullifiable_propositions.sort()

# Join the list into a single string with commas and no spaces
result = ",".join(nullifiable_propositions)

# Print the final result
print(result)