import math

def solve_spy_notebook():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """
    # Step 1: Determine the number of states per page
    num_fold_types = 3  # upper corner, lower corner, vertical half
    
    # States with 0 folds
    states_0_folds = 1
    
    # States with 1 fold (choose 1 type out of 3)
    states_1_fold = math.comb(num_fold_types, 1)
    
    # States with 2 folds (choose 2 different types out of 3)
    states_2_folds = math.comb(num_fold_types, 2)
    
    total_states_per_page = states_0_folds + states_1_fold + states_2_folds
    
    print(f"Each page can be in one of {total_states_per_page} states (1 for no fold + {states_1_fold} for one fold + {states_2_folds} for two folds).")
    
    # Step 2: Determine the number of unique observations
    num_sizes = 5
    num_times = 6
    total_observation_types = num_sizes * num_times
    
    print(f"There are {num_sizes} sizes and {num_times} time slots, making {total_observation_types} unique observation types to record.")
    print("-" * 20)

    # Step 3: Define notebook structure
    first_ordered_pages = 10
    last_ordered_pages = 10
    total_ordered_pages = first_ordered_pages + last_ordered_pages
    unordered_pages = 100 - total_ordered_pages

    # Step 4: Calculate observations for the 20 ordered pages
    print("For the first 10 and last 10 pages, James remembers the order.")
    pages_per_ordered_obs = 0
    p = 1
    while True:
        combinations = total_states_per_page ** p
        if combinations >= total_observation_types:
            pages_per_ordered_obs = p
            break
        p += 1
    
    print(f"To encode {total_observation_types} types, he needs an ordered sequence of {pages_per_ordered_obs} pages, giving {total_states_per_page}^{pages_per_ordered_obs}={combinations} possibilities.")
    
    num_ordered_obs = total_ordered_pages // pages_per_ordered_obs
    print(f"From the {total_ordered_pages} ordered pages, he can record {total_ordered_pages} / {pages_per_ordered_obs} = {num_ordered_obs} observations.")
    print("-" * 20)

    # Step 5: Calculate observations for the 80 unordered pages
    print(f"For the middle {unordered_pages} pages, the order is unknown.")
    pages_per_unordered_obs = 0
    p = 1
    while True:
        # Using combinations with repetition formula: C(n+k-1, k)
        # n = total_states_per_page, k = p (number of pages)
        combinations = math.comb(total_states_per_page + p - 1, p)
        if combinations >= total_observation_types:
            pages_per_unordered_obs = p
            break
        p += 1
        
    print(f"To encode {total_observation_types} types without order, he needs a set of {pages_per_unordered_obs} pages, giving C({total_states_per_page}+{pages_per_unordered_obs}-1, {pages_per_unordered_obs})={combinations} possibilities.")
    
    num_unordered_obs = unordered_pages // pages_per_unordered_obs
    print(f"From the {unordered_pages} unordered pages, he can record floor({unordered_pages} / {pages_per_unordered_obs}) = {num_unordered_obs} observations.")
    print("-" * 20)

    # Step 6: Sum the results
    total_observations = num_ordered_obs + num_unordered_obs
    print("The highest number of observations he can record is the sum from both parts.")
    print(f"Final Calculation: {num_ordered_obs} + {num_unordered_obs} = {total_observations}")

solve_spy_notebook()
<<<36>>>