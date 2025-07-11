import scipy.stats

# --- Define Test Parameters ---
# The initial significance level (alpha) for the test.
original_alpha = 0.05
# The number of statistical tests being performed.
num_tests = 5
# The observations obtained in one session for each group.
base_group1 = [1, 2, 3, 4, 5]
base_group2 = [3, 3, 3, 2, 2]

# --- Step 1: Calculate Bonferroni Corrected Alpha ---
# This is the new p-value threshold we need to beat.
corrected_alpha = original_alpha / num_tests

print("--- Goal: Find Minimum Sample Size for Mann-Whitney U Test ---")
print(f"Original p-value threshold: {original_alpha}")
print(f"Number of tests: {num_tests}")
print(f"Bonferroni corrected p-value threshold: {original_alpha} / {num_tests} = {corrected_alpha:.2f}\n")

# --- Step 2 & 3: Iterate and Test ---
# We start with 1 session and increase until significance is reached.
repetition_factor = 0
p_value = 1 # Initialize p-value to a value > alpha

while p_value >= corrected_alpha:
    repetition_factor += 1
    
    # Create the data groups by repeating the base observations
    current_group1 = base_group1 * repetition_factor
    current_group2 = base_group2 * repetition_factor

    # Perform the two-sided Mann-Whitney U test
    u_statistic, p_value = scipy.stats.mannwhitneyu(
        current_group1,
        current_group2,
        alternative='two-sided'
    )
    
    # Print the status of each iteration
    n_per_group = len(current_group1)
    print(f"Testing with {repetition_factor} session(s) ({n_per_group} observations/group)... p-value = {p_value:.5f}")

# --- Step 4: Report the Final Result ---
# The loop stops when p_value < corrected_alpha. The final values are from the last iteration.
final_n_per_group = len(current_group1)

print("\n--- Significance Achieved ---")
print(f"Minimum number of observations required per group: {final_n_per_group}")
print("\n--- Final Equation Check ---")
print("The condition for significance is: p-value < corrected alpha")
print(f"Our final result is: {p_value:.5f} < {corrected_alpha:.2f}")
print(f"This condition is met, so {final_n_per_group} is the minimum number of observations per group.")

<<<30>>>