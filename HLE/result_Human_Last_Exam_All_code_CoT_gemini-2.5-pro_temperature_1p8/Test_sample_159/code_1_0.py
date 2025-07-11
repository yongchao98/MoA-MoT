import math

# Define the constants from the problem description
NUM_PAGES_IN_NOTEBOOK = 100
NUM_SIZE_LEVELS = 5  # few, small, medium, large, huge
NUM_TIME_SLOTS = 6   # 12am, 4am, 8am, 12pm, 4pm, 8pm

print("Step 1: Calculate the number of distinct states for a single page.")
# The available folds are: Right Upper corner (U), Right Lower corner (L), Vertical half (V).
# At most two folds are allowed per page.

# Number of states with 0 folds:
# - No folds: 1 state
states_with_0_folds = 1

# Number of states with 1 fold:
# - {U}, {L}, or {V}: 3 distinct states
states_with_1_fold = 3

# Number of states with 2 folds:
# We choose 2 distinct folds from the 3 available types.
# - {U, L}, {U, V}, {L, V}: 3 distinct combinations
states_with_2_folds = 3

# Total states per page is the sum of states for 0, 1, and 2 folds.
total_states_per_page = states_with_0_folds + states_with_1_fold + states_with_2_folds
print(f"A page can have 0 folds (1 state), 1 fold (3 states), or 2 folds (3 states).")
print(f"Total states per page = {states_with_0_folds} + {states_with_1_fold} + {states_with_2_folds} = {total_states_per_page}\n")

print("Step 2: Calculate the number of distinct observation types.")
# Each observation is a combination of a size and a time.
total_observation_types = NUM_SIZE_LEVELS * NUM_TIME_SLOTS
print(f"There are {NUM_SIZE_LEVELS} size levels and {NUM_TIME_SLOTS} time slots.")
print(f"Total unique observation types = {NUM_SIZE_LEVELS} * {NUM_TIME_SLOTS} = {total_observation_types}\n")

print("Step 3: Determine the pages needed to encode one observation.")
# We need to find the smallest integer 'k' (pages) such that the number of states
# we can represent with 'k' pages is enough to cover all observation types.
# The condition is: (total_states_per_page)^k >= total_observation_types
pages_per_observation = math.ceil(math.log(total_observation_types) / math.log(total_states_per_page))

print(f"We need to find the smallest integer 'k' where {total_states_per_page}^k >= {total_observation_types}.")
# Showing the check for k=1 and k=2
print(f"For k=1, {total_states_per_page}^1 = {total_states_per_page}, which is less than {total_observation_types}.")
print(f"For k=2, {total_states_per_page}^2 = {total_states_per_page**2}, which is greater than or equal to {total_observation_types}.")
print(f"Therefore, {pages_per_observation} pages are needed to record a single observation.\n")

print("Step 4: Calculate the highest number of observations James can record.")
# The total number of observations is the notebook's total pages divided by the pages needed per observation.
# We use integer division as we can only record a whole number of observations.
max_observations = NUM_PAGES_IN_NOTEBOOK // pages_per_observation
print(f"With a {NUM_PAGES_IN_NOTEBOOK}-page notebook and {pages_per_observation} pages required per observation, the total is calculated as follows:")
print(f"Highest number of observations = Total Pages / Pages per Observation")
print(f"                               = {NUM_PAGES_IN_NOTEBOOK} / {pages_per_observation} = {max_observations}")
<<<50>>>