import math

def solve_spy_notebook_problem():
    """
    Calculates the highest number of observations James can record in his notebook.
    """

    # Step 1: Define the constants from the problem.
    order_matters_pages = 20  # First 10 and last 10 pages
    order_doesnt_matter_pages = 80  # Middle 80 pages
    observation_outcomes = 5  # 5 size levels: few, small, medium, large, huge

    # Step 2: Calculate the number of states for the "no order" pages.
    # For these pages, only the final physical configuration matters.
    # States are defined by the set of unique folds applied.
    # 0 folds: 1 state (unchanged)
    # 1 fold: 3 states {U}, {L}, {V}
    # 2 folds: 4 distinct states based on combinations: {U,L}, {U,V}, {L,V}, and {V,V} (quartered page)
    # The states {U,U} and {L,L} are physically identical to {U} and {L} respectively.
    no_order_states = 1 + 3 + 4

    # Step 3: Calculate the number of states for the "order matters" pages.
    # For these pages, the sequence of operations is the information.
    # 0 folds: 1 state (empty sequence)
    # 1 fold: 3 states (U), (L), (V)
    # 2 folds: 3 * 3 = 9 states. All sequences ((U,U), (U,L), etc.) are distinguishable
    # because James remembers the sequence of operations he performed.
    order_states = 1 + 3 + (3 * 3)

    # Step 4: Formulate the inequality to find the number of observations N.
    # The total number of ways to record N observations is observation_outcomes^N.
    # The total information capacity of the notebook is order_states^order_matters_pages * no_order_states^order_doesnt_matter_pages.
    # We need to find the largest integer N such that:
    # observation_outcomes^N <= order_states^order_matters_pages * no_order_states^order_doesnt_matter_pages
    #
    # Taking the logarithm of both sides:
    # N * log(observation_outcomes) <= order_matters_pages * log(order_states) + order_doesnt_matter_pages * log(no_order_states)
    #
    # N <= (order_matters_pages * log(order_states) + order_doesnt_matter_pages * log(no_order_states)) / log(observation_outcomes)

    # Step 5: Calculate the values and solve for N.
    log_order_states = math.log(order_states)
    log_no_order_states = math.log(no_order_states)
    log_observation_outcomes = math.log(observation_outcomes)

    numerator = (order_matters_pages * log_order_states) + (order_doesnt_matter_pages * log_no_order_states)
    denominator = log_observation_outcomes
    
    max_n = numerator / denominator
    result = math.floor(max_n)
    
    # Step 6: Print the thinking process and the final answer.
    print(f"The calculation is based on the inequality: {observation_outcomes}^N <= {order_states}^{order_matters_pages} * {no_order_states}^{order_doesnt_matter_pages}")
    print(f"This is solved for N by taking the logarithm: N <= ({order_matters_pages} * log({order_states}) + {order_doesnt_matter_pages} * log({no_order_states})) / log({observation_outcomes})")
    print("\n--- Values for the Equation ---")
    print(f"Number of pages where folding order matters: {order_matters_pages}")
    print(f"States per 'order matters' page: {order_states}")
    print(f"Number of pages where folding order does not matter: {order_doesnt_matter_pages}")
    print(f"States per 'no order' page: {no_order_states}")
    print(f"Possible outcomes for each observation: {observation_outcomes}")
    print("---------------------------------")
    print(f"\nThe highest number of observations James can record is: {result}")

solve_spy_notebook_problem()