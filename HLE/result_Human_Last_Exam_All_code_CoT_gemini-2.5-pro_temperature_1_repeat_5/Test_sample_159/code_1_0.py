import math

def calculate_max_observations():
    """
    Calculates the maximum number of observations James can record in his notebook.
    """
    # Step 1: Define problem parameters and calculate unique observation types
    num_sizes = 5  # few, small, medium, large, huge
    num_times = 6   # 12am, 4am, 8am, 12pm, 4pm, 8pm
    observation_types_needed = num_sizes * num_times
    
    total_pages = 100
    num_ordered_pages = 10 + 10
    num_unordered_pages = total_pages - num_ordered_pages

    print("--- Step 1: Analyzing the Observation ---")
    print(f"An observation consists of 1 of {num_sizes} sizes and 1 of {num_times} times.")
    print(f"Total unique observation types to encode: {num_sizes} * {num_times} = {observation_types_needed}\n")

    # Step 2: Calculate the number of states per page
    num_base_folds = 3  # Right Upper, Right Lower, Vertical

    # Unordered Pages (order doesn't matter, so we use combinations)
    # 0 folds: C(3,0) = 1 state (no fold)
    # 1 fold: C(3,1) = 3 states
    # 2 folds: C(3,2) = 3 states
    unordered_page_states = 1 + math.comb(num_base_folds, 1) + math.comb(num_base_folds, 2)
    
    # Ordered Pages (order matters, so we use permutations)
    # 0 folds: 1 state
    # 1 fold: P(3,1) = 3 states
    # 2 folds: P(3,2) = 6 states
    ordered_page_states = 1 + math.perm(num_base_folds, 1) + math.perm(num_base_folds, 2)
    
    print("--- Step 2: Calculating Page States ---")
    print(f"For the {num_unordered_pages} middle pages, fold order doesn't matter.")
    print(f"States for an unordered page = (0 folds) + (1 fold) + (2 folds) = 1 + {math.comb(num_base_folds, 1)} + {math.comb(num_base_folds, 2)} = {unordered_page_states}")
    print(f"For the {num_ordered_pages} first/last pages, fold order matters.")
    print(f"States for an ordered page = (0 folds) + (1 fold) + (2 folds) = 1 + {math.perm(num_base_folds, 1)} + {math.perm(num_base_folds, 2)} = {ordered_page_states}\n")

    # Step 3: Determine how to combine pages for one observation
    # We need a combination of pages that can represent at least 30 states.
    # Using 1 page is not enough (7 < 30 and 10 < 30). We need 2 pages.
    print("--- Step 3: Pages per Observation ---")
    print(f"To encode {observation_types_needed} states, we need to combine pages:")
    print(f" - 2 unordered pages: {unordered_page_states} * {unordered_page_states} = {unordered_page_states**2} states (Sufficient)")
    print(f" - 1 ordered & 1 unordered page: {ordered_page_states} * {unordered_page_states} = {ordered_page_states * unordered_page_states} states (Sufficient)")
    print(f" - 2 ordered pages: {ordered_page_states} * {ordered_page_states} = {ordered_page_states**2} states (Sufficient)\n")

    # Step 4: Solve the resource allocation problem to maximize observations
    # Let x = # observations using 2 unordered pages
    # Let y = # observations using 1 ordered and 1 unordered page
    # Let z = # observations using 2 ordered pages
    # Maximize N = x + y + z
    # Subject to constraints:
    # 2*x + y <= 80 (unordered pages)
    # y + 2*z <= 20 (ordered pages)
    # By solving this linear system, we find the maximum N is 50.
    # One solution is y=0, which gives x<=40 and z<=10.
    # To maximize N = x+z, we choose x=40 and z=10.
    x = 40
    y = 0
    z = 10
    max_observations = x + y + z

    print("--- Step 4: Maximizing Total Observations ---")
    print("We have a limited supply of pages. We can set up an equation to maximize the number of observations.")
    print("Let 'x' be observations using 2 unordered pages, 'y' using one of each, and 'z' using 2 ordered pages.")
    print(f"The optimal allocation is:")
    print(f"  - {x} observations using 2 unordered pages each (uses {x*2} unordered pages).")
    print(f"  - {y} observations using 1 of each type (uses {y} unordered and {y} ordered pages).")
    print(f"  - {z} observations using 2 ordered pages each (uses {z*2} ordered pages).")
    print("\nThis allocation uses all available pages and maximizes the result.")
    print(f"Final Calculation: {x} + {y} + {z} = {max_observations}")
    print("\n-------------------------------------------------")
    print(f"The highest number of observations is: {max_observations}")
    print("-------------------------------------------------")


if __name__ == "__main__":
    calculate_max_observations()