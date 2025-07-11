def solve_spy_notebook():
    """
    Calculates the maximum number of observations James can record.
    """

    # 1. Deconstruct the Observation
    # An observation has a size and a time.
    num_sizes = 5  # few, small, medium, large, huge
    num_times = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm

    # 2. Analyze the Recording Medium
    total_pages = 100
    
    # 3. Calculate Page States
    # A page can be unchanged (1 state), have one fold (3 types), or two folds.
    # The minimum number of states a page can have is when the folding order doesn't matter.
    # This applies to the middle 80 pages.
    # - 0 folds: 1 state
    # - 1 fold: 3 states (upper corner, lower corner, vertical)
    # - 2 folds (order doesn't matter): Choosing 2 distinct folds from 3 types is C(3,2) = 3 states.
    # Total states for such a page = 1 + 3 + 3 = 7.
    min_states_per_page = 7

    # 4. Devise an Encoding Strategy
    # To record one observation, we need to record one size (5 options) and one time (6 options).
    # We can use one page for size and one page for time.
    # A page needs at least 5 states to encode size, and 6 states to encode time.
    # Since even the pages with the minimum states (7) have enough capacity (7 > 5 and 7 > 6),
    # any page can be used for either task.
    # Therefore, each observation requires two pages.
    pages_per_observation = 2

    # 5. Calculate the Total Capacity
    # The maximum number of observations is the total pages divided by the pages needed per observation.
    max_observations = total_pages // pages_per_observation

    print("Step 1: Determine the information for one observation.")
    print(f"An observation requires encoding a size ({num_sizes} options) and a time ({num_times} options).")
    print("\nStep 2: Determine the pages needed per observation.")
    print(f"A single page has at least {min_states_per_page} states, which is enough to encode either size or time.")
    print(f"Therefore, one observation requires {pages_per_observation} pages (one for size, one for time).")
    print("\nStep 3: Calculate the maximum number of observations for the whole notebook.")
    print("The final equation is:")
    print(f"{total_pages} (total pages) / {pages_per_observation} (pages per observation) = {max_observations}")

solve_spy_notebook()