import math

def solve_spy_notebook_problem():
    """
    Calculates the maximum number of observations James can record in his notebook.

    This function follows these steps:
    1.  Calculates the number of states for a single observation (Size x Time).
    2.  Calculates the number of states for the two types of pages:
        a. 'Unordered' pages where only the final fold combination matters.
        b. 'Ordered' pages where the sequence of folds matters.
    3.  Calculates the total information capacity of the notebook in bits (log base 2),
        though any log base works for the intermediate calculation.
    4.  Uses the information capacity to find the maximum number of observations (N)
        that can be encoded by solving the inequality:
        Total_Notebook_States >= Total_Observation_States^N
    """
    # Step 1: Define the parameters for an observation
    num_sizes = 5  # {few, small, medium, large, huge}
    num_times = 6  # {12am, 4am, 8am, 12pm, 4pm, 8pm}
    states_per_observation = num_sizes * num_times

    # Step 2: Define the parameters for the notebook and folds
    total_pages = 100
    num_ordered_pages = 10 + 10  # First 10 and last 10
    num_unordered_pages = total_pages - num_ordered_pages
    num_fold_types = 3  # {right upper corner, right lower corner, vertical half}

    # Calculate states for an UNORDERED page (order of 2 folds doesn't matter)
    # 0 folds: 1 state (unchanged)
    # 1 fold: C(3, 1) = 3 states
    # 2 folds: This is combinations with repetition. C(n+k-1, k) for n=3, k=2.
    # C(3+2-1, 2) = C(4, 2) = (4*3)/(2*1) = 6 states.
    states_per_unordered_page = 1 + 3 + 6

    # Calculate states for an ORDERED page (order of 2 folds matters)
    # 0 folds: 1 state (unchanged)
    # 1 fold: P(3, 1) = 3 states
    # 2 folds: This is permutations with repetition. n^k for n=3, k=2.
    # 3^2 = 9 states.
    states_per_ordered_page = 1 + 3 + 9

    # Step 3 & 4: Use logarithms to solve for N in the inequality:
    # (states_ordered^num_ordered) * (states_unordered^num_unordered) >= states_observation^N
    #
    # Taking the log of both sides:
    # num_ordered*log(states_ordered) + num_unordered*log(states_unordered) >= N*log(states_observation)
    #
    # N <= (num_ordered*log(states_ordered) + num_unordered*log(states_unordered)) / log(states_observation)

    log_total_notebook_capacity = (num_ordered_pages * math.log(states_per_ordered_page) +
                                   num_unordered_pages * math.log(states_per_unordered_page))

    log_states_per_observation = math.log(states_per_observation)

    max_n_float = log_total_notebook_capacity / log_states_per_observation
    max_observations = math.floor(max_n_float)

    # Print the explanation and the final equation
    print("### Step-by-Step Calculation ###\n")
    print(f"1. States per page (where fold order is remembered):")
    print(f"   1 (no fold) + 3 (one fold) + 9 (two ordered folds) = {states_per_ordered_page} states")
    print(f"   There are {num_ordered_pages} such pages.\n")

    print(f"2. States per page (where fold order is forgotten):")
    print(f"   1 (no fold) + 3 (one fold) + 6 (two unordered folds) = {states_per_unordered_page} states")
    print(f"   There are {num_unordered_pages} such pages.\n")
    
    print(f"3. States per observation:")
    print(f"   {num_sizes} (sizes) * {num_times} (times) = {states_per_observation} states\n")

    print(f"4. Finding the maximum number of observations (N):")
    print(f"   The problem is to solve: ({states_per_ordered_page}^{num_ordered_pages}) * ({states_per_unordered_page}^{num_unordered_pages}) >= {states_per_observation}^N")
    print(f"   Using logarithms, we solve for N:")
    print(f"   N <= ({num_ordered_pages} * log({states_per_ordered_page}) + {num_unordered_pages} * log({states_per_unordered_page})) / log({states_per_observation})")
    print(f"   N <= {max_n_float:.2f}\n")
    
    print(f"Since N must be a whole number, the highest number of observations James can record is the integer part.")
    print(f"Highest number of observations = {max_observations}")
    
    # Final answer in the required format
    print(f"\n<<<{max_observations}>>>")

solve_spy_notebook_problem()