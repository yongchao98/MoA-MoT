import math

def solve_spy_notebook():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """
    # --- Constants based on the problem statement ---
    TOTAL_PAGES = 100
    NUM_SIZES = 5  # few, small, medium, large, huge
    NUM_TIMES = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
    NUM_BASIC_FOLDS = 3 # Upper corner, Lower corner, Vertical half

    print("Step 1: Calculate the number of states a single page can represent.")
    
    # For the first 10 and last 10 pages, the order of two folds matters.
    # States = 1 (no fold) + 3 (one fold) + Permutations of 2 from 3 (P(3,2))
    # P(3,2) = 3! / (3-2)! = 6
    no_fold_states = 1
    one_fold_states = NUM_BASIC_FOLDS
    two_folds_order_matters = math.factorial(NUM_BASIC_FOLDS) // math.factorial(NUM_BASIC_FOLDS - 2)
    states_order_matters = no_fold_states + one_fold_states + two_folds_order_matters
    print(f"For pages where fold order matters (20 pages total):")
    print(f"{no_fold_states} (no fold) + {one_fold_states} (one fold) + {two_folds_order_matters} (two folds) = {states_order_matters} states.")

    # For the middle 80 pages, the order of two folds does not matter.
    # States = 1 (no fold) + 3 (one fold) + Combinations of 2 from 3 (C(3,2))
    # C(3,2) = 3! / (2! * (3-2)!) = 3
    two_folds_order_doesnt_matter = math.factorial(NUM_BASIC_FOLDS) // (math.factorial(2) * math.factorial(NUM_BASIC_FOLDS - 2))
    states_order_doesnt_matter = no_fold_states + one_fold_states + two_folds_order_doesnt_matter
    print(f"For pages where fold order does not matter (80 pages total):")
    print(f"{no_fold_states} (no fold) + {one_fold_states} (one fold) + {two_folds_order_doesnt_matter} (two folds) = {states_order_doesnt_matter} states.\n")

    print("Step 2: Calculate the number of different observations to be recorded.")
    observation_types = NUM_SIZES * NUM_TIMES
    print(f"Total unique observation types = {NUM_SIZES} (sizes) * {NUM_TIMES} (times) = {observation_types} types.\n")

    print("Step 3: Determine the minimum number of pages needed per observation.")
    # A single page can represent at most 10 states, which is not enough for 30 observation types.
    print(f"One page can represent at most {states_order_matters} states, which is less than the required {observation_types} states.")
    
    # We check if two pages are sufficient. We use the minimum number of states per page
    # to ensure any pair of pages will work.
    min_states_per_page = states_order_doesnt_matter
    pages_per_observation = 2
    states_with_two_pages = min_states_per_page ** pages_per_observation
    print(f"Using {pages_per_observation} pages, we can represent at least {min_states_per_page} * {min_states_per_page} = {states_with_two_pages} states.")
    print(f"Since {states_with_two_pages} > {observation_types}, two pages are sufficient to record any observation.\n")

    print("Step 4: Calculate the maximum number of observations possible.")
    max_observations = TOTAL_PAGES // pages_per_observation
    print(f"With {TOTAL_PAGES} total pages and {pages_per_observation} pages per observation, James can record a maximum of:")
    print(f"{TOTAL_PAGES} / {pages_per_observation} = {max_observations} observations.")

solve_spy_notebook()
<<<50>>>