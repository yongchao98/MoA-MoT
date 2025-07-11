import math

# Step 1: Calculate the number of states for a page without folding order memory.
# There are 3 distinct types of single folds.
# - 0 folds: 1 state (unchanged).
# - 1 fold: Choosing 1 fold type from 3 is C(3, 1) = 3 states.
# - 2 folds: Choosing 2 distinct fold types from 3 is C(3, 2) = 3 states.
states_unordered = 1 + 3 + 3
num_unordered_pages = 80

# Step 2: Calculate the number of states for a page with folding order memory.
# - 0 folds: 1 state (unchanged).
# - 1 fold: The sequence has one fold, so there are 3 states.
# - 2 folds: The sequence has two folds. Each can be any of the 3 types, so 3 * 3 = 9 states.
states_ordered = 1 + 3 + 9
num_ordered_pages = 20

# Step 3: Calculate the number of states needed for one observation.
# There are 5 possible sizes and 6 possible times.
num_sizes = 5
num_times = 6
states_per_observation = num_sizes * num_times

# Step 4: Formulate the problem to find the maximum number of observations (N).
# The total number of unique notebook configurations must be at least the
# total number of possible sequences of N observations.
# Inequality: states_per_observation^N <= states_ordered^num_ordered_pages * states_unordered^num_unordered_pages

# Step 5: Solve for N using logarithms.
# N * log(states_per_observation) <= num_ordered_pages * log(states_ordered) + num_unordered_pages * log(states_unordered)
# N <= (num_ordered_pages * log(states_ordered) + num_unordered_pages * log(states_unordered)) / log(states_per_observation)

# Perform the calculation
log_total_states = num_ordered_pages * math.log(states_ordered) + num_unordered_pages * math.log(states_unordered)
log_obs_states = math.log(states_per_observation)
max_N = log_total_states / log_obs_states
result = math.floor(max_N)

# Print the results step-by-step
print("Step 1: Determine the information capacity of each page type.")
print(f"A page without order memory can represent {states_unordered} states.")
print(f"A page with order memory can represent {states_ordered} states.")
print("\nStep 2: Determine the information required per observation.")
print(f"Each observation requires encoding one of {states_per_observation} possibilities ({num_sizes} sizes * {num_times} times).")
print("\nStep 3: Set up the equation to find the maximum number of observations (N).")
print("The total states of the notebook must be greater than or equal to the total states for N observations.")
print(f"This gives the inequality: {states_per_observation}^N <= {states_ordered}^{num_ordered_pages} * {states_unordered}^{num_unordered_pages}")
print("\nStep 4: Solve the equation using logarithms.")
print(f"The equation to solve for N is:")
print(f"N <= ({num_ordered_pages} * log({states_ordered}) + {num_unordered_pages} * log({states_unordered})) / log({states_per_observation})")
print("\nStep 5: Calculate the final result.")
print(f"N <= {max_N:.4f}")
print(f"\nThe highest number of observations James can record is the largest integer less than or equal to {max_N:.4f}, which is {result}.")
