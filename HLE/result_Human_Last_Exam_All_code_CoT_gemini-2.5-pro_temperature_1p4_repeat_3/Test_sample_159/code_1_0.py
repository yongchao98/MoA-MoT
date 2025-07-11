import math

def solve_spy_notebook_problem():
    """
    Calculates the highest number of observations a spy can record in a notebook.
    The calculation is based on the principles of combinatorics and information theory.
    """

    # --- Step 1: Define the fundamental quantities ---

    # A page can be in one of the following states:
    # 1. Unchanged (no folds)
    # 2. One fold (upper corner, lower corner, or vertical)
    # 3. Two folds (combinations of the three single fold types: UL, UV, LV)
    states_per_page = 1 + 3 + 3  # 1 (no fold) + 3 (one fold) + 3 (two distinct folds)

    # An observation is a combination of a size estimate and a time.
    size_categories = 5  # few, small, medium, large, huge
    time_slots = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
    observation_types = size_categories * time_slots

    # The notebook's page structure.
    total_pages = 100
    ordered_pages_per_end = 10
    num_ordered_pages = 2 * ordered_pages_per_end
    num_unordered_pages = total_pages - num_ordered_pages

    print("### Plan Summary ###")
    print("1. Determine the number of states per page and types of observations.")
    print("2. Calculate the total information capacity of the notebook (Total States).")
    print("3. Solve for N in (Observation Types)^N <= (Total States) using logarithms.")
    print("-" * 30)

    # --- Step 2: Calculate the notebook's total information capacity ---

    # For the first 10 and last 10 pages, the order is remembered.
    # The number of configurations is a permutation with repetition: states^pages
    # We use logarithms to handle the large numbers.
    log_capacity_ordered = num_ordered_pages * math.log(states_per_page)

    # For the middle 80 pages, order is not remembered.
    # This is a combination with repetition problem (stars and bars).
    # n = num_unordered_pages (items), k = states_per_page (bins)
    # Number of combinations = C(n + k - 1, k - 1)
    capacity_unordered = math.comb(num_unordered_pages + states_per_page - 1, states_per_page - 1)
    log_capacity_unordered = math.log(capacity_unordered)

    # The total capacity is the sum of the log capacities.
    log_total_capacity = log_capacity_ordered + log_capacity_unordered

    # --- Step 3: Find the maximum number of observations (N) ---

    # We need to solve for the largest integer N in:
    # observation_types^N <= total_capacity
    # N * log(observation_types) <= log_total_capacity
    # N <= log_total_capacity / log(observation_types)
    log_observation_types = math.log(observation_types)
    max_n = log_total_capacity / log_observation_types
    final_answer = math.floor(max_n)

    # --- Step 4: Display the results with the final equation ---

    print("### Step-by-Step Calculation ###\n")
    print(f"Number of states per page: {states_per_page}")
    print(f"Number of observation types ({size_categories} sizes * {time_slots} times): {observation_types}\n")
    
    print(f"The notebook has {num_ordered_pages} pages where order matters and {num_unordered_pages} where it does not.")
    print(f"Capacity of ordered pages = {states_per_page}^{num_ordered_pages}")
    print(f"Capacity of unordered pages = C({num_unordered_pages} + {states_per_page} - 1, {states_per_page} - 1) = C(86, 6) = {capacity_unordered}\n")

    print("### Final Equation ###")
    print(f"We solve for N: {observation_types}^N <= {states_per_page}^{num_ordered_pages} * {capacity_unordered}")
    print("Using logarithms:")
    print(f"N <= (log({states_per_page}^{num_ordered_pages}) + log({capacity_unordered})) / log({observation_types})")
    print(f"N <= ({num_ordered_pages} * log({states_per_page}) + log({capacity_unordered})) / log({observation_types})")
    print(f"N <= ({num_ordered_pages} * {math.log(states_per_page):.4f} + {log_capacity_unordered:.4f}) / {log_observation_types:.4f}")
    print(f"N <= ({log_capacity_ordered:.4f} + {log_capacity_unordered:.4f}) / {log_observation_types:.4f}")
    print(f"N <= {log_total_capacity:.4f} / {log_observation_types:.4f}")
    print(f"N <= {max_n:.4f}\n")
    
    print("The highest number of observations is the integer part of this result.")
    print(f"Highest number of observations: {final_answer}")
    
# Run the solver
solve_spy_notebook_problem()