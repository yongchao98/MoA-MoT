import math

def solve_spy_notebook():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Calculate the total number of unique observation types.
    num_sizes = 5  # few, small, medium, large, huge
    num_times = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
    total_observation_types = num_sizes * num_times

    print("Step 1: Calculate the number of unique observations to record.")
    print(f"An observation consists of a size ({num_sizes} options) and a time ({num_times} options).")
    print(f"Total unique observation types = {num_sizes} * {num_times} = {total_observation_types}\n")

    # Step 2: Calculate the information capacity (number of states) of a single page.
    # There are 3 distinct fold types: Upper Corner, Lower Corner, Vertical.
    num_fold_types = 3

    # For the 80 middle pages, the order of folds does not matter (Combinations).
    # States = 1 (for 0 folds) + C(3,1) (for 1 fold) + C(3,2) (for 2 folds)
    states_normal_page = 1 + math.comb(num_fold_types, 1) + math.comb(num_fold_types, 2)
    
    print("Step 2: Calculate the number of states a single page can represent.")
    print("For the 80 middle pages, folding order doesn't matter.")
    print(f"Number of states for a normal page = 1 (no fold) + {math.comb(num_fold_types, 1)} (one fold) + {math.comb(num_fold_types, 2)} (two folds) = {states_normal_page}\n")

    # For the 20 special pages, the order of folds matters (Permutations).
    # States = 1 (for 0 folds) + P(3,1) (for 1 fold) + P(3,2) (for 2 folds)
    states_special_page = 1 + math.perm(num_fold_types, 1) + math.perm(num_fold_types, 2)

    print("For the 20 special pages (first 10, last 10), folding order matters.")
    print(f"Number of states for a special page = 1 (no fold) + {math.perm(num_fold_types, 1)} (one fold) + {math.perm(num_fold_types, 2)} (two folds) = {states_special_page}\n")

    # Step 3: Determine the minimum number of pages required per observation.
    print("Step 3: Determine the minimum number of pages per observation.")
    max_states_one_page = max(states_normal_page, states_special_page)
    print(f"A single page can represent at most {max_states_one_page} states.")
    print(f"Since {max_states_one_page} < {total_observation_types}, one page is not enough.\n")

    # Check if two pages are sufficient. We check the worst case (two normal pages).
    states_two_pages = states_normal_page * states_normal_page
    pages_per_observation = 2
    print(f"Checking with {pages_per_observation} pages:")
    print(f"Two normal pages can represent {states_normal_page} * {states_normal_page} = {states_two_pages} states.")
    print(f"Since {states_two_pages} > {total_observation_types}, two pages are sufficient to record any observation.\n")

    # Step 4: Calculate the maximum number of observations for the entire notebook.
    total_pages = 100
    max_observations = total_pages // pages_per_observation

    print("Step 4: Calculate the highest number of observations possible.")
    print(f"The notebook has a total of {total_pages} pages.")
    print(f"Each observation requires {pages_per_observation} pages.")
    print("The final calculation is:")
    print(f"Total Pages / Pages per Observation = {total_pages} / {pages_per_observation} = {max_observations}")

solve_spy_notebook()