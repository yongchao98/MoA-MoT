import math

# Step 1: Define the problem parameters
num_sizes = 5
num_times = 6
total_pages = 100
ordered_pages_count = 10 + 10
unordered_pages_count = total_pages - ordered_pages_count

# Calculate the number of unique states required for one observation
states_per_observation = num_sizes * num_times

# Step 2: Calculate the information capacity (number of states) per page type
# For 'ordered' pages, order matters (Permutations)
# 1 (no fold) + P(3,1) + P(3,2)
states_per_ordered_page = 1 + 3 + (math.factorial(3) // math.factorial(3 - 2))

# For 'unordered' pages, order does not matter (Combinations)
# 1 (no fold) + C(3,1) + C(3,2)
states_per_unordered_page = 1 + 3 + (math.factorial(3) // (math.factorial(2) * math.factorial(3 - 2)))

# Step 3: Calculate the number of pages required to encode one observation
# We need to find k such that states^k >= states_per_observation
# For ordered pages: 10^k >= 30 -> k=2
# For unordered pages: 7^k >= 30 -> k=2
# We can calculate this systematically, but for this problem, it's clearly 2.
pages_per_obs_ordered = 2
pages_per_obs_unordered = 2

# Step 4: Calculate the total number of observations
num_obs_from_ordered = ordered_pages_count // pages_per_obs_ordered
num_obs_from_unordered = unordered_pages_count // pages_per_obs_unordered
total_observations = num_obs_from_ordered + num_obs_from_unordered

# Print the final equation with all the numbers
print("The highest number of observations is calculated as follows:")
print(f"Total Observations = (Total Ordered Pages / Pages per Observation) + (Total Unordered Pages / Pages per Observation)")
print(f"Total Observations = ({ordered_pages_count} / {pages_per_obs_ordered}) + ({unordered_pages_count} / {pages_per_obs_unordered})")
print(f"Total Observations = {num_obs_from_ordered} + {num_obs_from_unordered}")
print(f"Total Observations = {total_observations}")
