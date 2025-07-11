import math

# Step 1: Calculate the number of states per page type.

# There are 3 distinct types of single folds James can make.
num_fold_types = 3

# For the 80 middle pages, the order of folds does not matter (combinations).
# 0 folds: 1 state (unchanged)
states_0_folds_middle = 1
# 1 fold: C(3, 1)
states_1_fold_middle = math.comb(num_fold_types, 1)
# 2 folds: C(3, 2)
states_2_folds_middle = math.comb(num_fold_types, 2)
total_states_middle_page = states_0_folds_middle + states_1_fold_middle + states_2_folds_middle

print(f"A single middle page (where fold order doesn't matter) can have {total_states_middle_page} distinct states.")

# For the 20 special pages (first 10, last 10), the order of folds matters (permutations).
# 0 folds: 1 state (unchanged)
states_0_folds_special = 1
# 1 fold: P(3, 1)
states_1_fold_special = math.perm(num_fold_types, 1)
# 2 folds: P(3, 2)
states_2_folds_special = math.perm(num_fold_types, 2)
total_states_special_page = states_0_folds_special + states_1_fold_special + states_2_folds_special

print(f"A single special page (where fold order matters) can have {total_states_special_page} distinct states.")
print("-" * 20)

# Step 2: Calculate the number of states required per observation.
num_size_options = 5  # few, small, medium, large, huge
num_time_options = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
states_per_observation = num_size_options * num_time_options

print(f"Number of size options: {num_size_options}")
print(f"Number of time options: {num_time_options}")
print(f"Total unique observations to record: {num_size_options} * {num_time_options} = {states_per_observation} states.")
print("-" * 20)

# Step 3: Determine the minimum number of pages per observation.
print(f"One page is not enough because the maximum states on a single page is {total_states_special_page}, which is less than {states_per_observation}.")

# Let's check if 2 pages are enough. We check the weakest combination (two middle pages).
states_with_2_middle_pages = total_states_middle_page * total_states_middle_page
print(f"Using two middle pages, we can represent {total_states_middle_page} * {total_states_middle_page} = {states_with_2_middle_pages} states.")
print(f"Since {states_with_2_middle_pages} is greater than {states_per_observation}, 2 pages are sufficient to record one observation.")

pages_per_observation = 2
print(f"Minimum pages needed per observation: {pages_per_observation}")
print("-" * 20)

# Step 4: Calculate the maximum number of observations.
total_pages = 100
max_observations = total_pages // pages_per_observation

print("The highest number of observations is calculated as follows:")
print(f"Total Pages / Pages per Observation = {total_pages} / {pages_per_observation} = {max_observations}")
print("-" * 20)
print(f"Therefore, the highest number of observations James can record is {max_observations}.")
