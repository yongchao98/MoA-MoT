import math

def solve_spy_notebook():
    """
    Calculates the highest number of observations James can record.
    """
    # Step 1: Calculate the number of states per page.
    # A page can have 0, 1, or 2 folds.
    # The 3 types of folds are: right upper corner, right lower corner, vertical half.
    # - 0 folds: 1 state (unchanged).
    # - 1 fold: 3 states (one for each fold type).
    # - 2 folds: Combinations of choosing 2 distinct fold types from 3. C(3, 2) = 3.
    states_per_page = 1 + 3 + math.comb(3, 2)

    # Step 2: Calculate the number of states per observation.
    # 5 size options * 6 time options.
    size_options = 5
    time_options = 6
    states_per_observation = size_options * time_options

    # Step 3: Calculate the total information capacity of the notebook.
    # The notebook has two types of pages.
    # - 20 "ordered" pages (first 10 + last 10).
    # - 80 "unordered" pages (the middle ones).
    num_ordered_pages = 20
    num_unordered_pages = 80

    # For the ordered part, the total states are states_per_page ^ num_ordered_pages.
    # We use logarithms to handle the very large numbers.
    log_total_states = num_ordered_pages * math.log(states_per_page)

    # For the unordered part, we have a combinations with repetition problem.
    # It's equivalent to distributing n identical items (pages) into k distinct bins (states).
    # The formula is C(n + k - 1, k - 1).
    # Here, n = num_unordered_pages and k = states_per_page.
    unordered_combinations = math.comb(num_unordered_pages + states_per_page - 1, states_per_page - 1)
    log_total_states += math.log(unordered_combinations)

    # Step 4: Determine the maximum number of observations (N).
    # The inequality is: states_per_observation ^ N <= total_notebook_states.
    # Taking the log of both sides: N * log(states_per_observation) <= log_total_states.
    # N <= log_total_states / log(states_per_observation).
    log_states_per_obs = math.log(states_per_observation)
    max_n_float = log_total_states / log_states_per_obs
    
    # The number of observations must be an integer.
    max_n_integer = math.floor(max_n_float)

    # Output the components of the final equation as requested.
    print(f"This problem is solved by comparing the information capacity of the notebook with the information required for the observations.")
    print(f"The equation is: (States per Observation)^N <= (Total Notebook States)")
    print(f"\n1. States per Page: {states_per_page}")
    print(f"2. States per Observation: {states_per_observation}")
    print(f"3. Total Notebook States are calculated from:")
    print(f"   - {num_ordered_pages} ordered pages: {states_per_page}^{num_ordered_pages}")
    print(f"   - {num_unordered_pages} unordered pages: C({num_unordered_pages} + {states_per_page} - 1, {states_per_page} - 1) = {unordered_combinations}")
    print(f"\nSolving for N in the inequality:")
    print(f"{states_per_observation}^N <= {states_per_page}^{num_ordered_pages} * {unordered_combinations}")
    print(f"\nThe highest number of observations (N) James can record is: {max_n_integer}")


solve_spy_notebook()