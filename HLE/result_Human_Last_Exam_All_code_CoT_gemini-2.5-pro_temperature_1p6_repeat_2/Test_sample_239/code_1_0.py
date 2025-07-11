import math

# Define the input values for the Theatre Square problem.
# These values correspond to the specific case asked in question 4.
n = 4 * 10**9
m = 4 * 10**9
a = 1

# To find the number of flagstones needed along one dimension, we need to
# calculate the ceiling of the division (e.g., ceil(n/a)).
# A common and branchless way to do this with integer arithmetic is (x + y - 1) // y.

# Calculate the number of flagstones needed for the n-dimension.
num_stones_n = (n + a - 1) // a

# Calculate the number of flagstones needed for the m-dimension.
num_stones_m = (m + a - 1) // a

# The total number of flagstones is the product of the two dimensions.
total_stones = num_stones_n * num_stones_m

# Print the final equation with all its numbers as requested.
print(f"To cover the n-dimension, we need: ({n} + {a} - 1) // {a} = {num_stones_n} flagstones.")
print(f"To cover the m-dimension, we need: ({m} + {a} - 1) // {a} = {num_stones_m} flagstones.")
print(f"The total number of flagstones is: {num_stones_n} * {num_stones_m} = {total_stones}")
