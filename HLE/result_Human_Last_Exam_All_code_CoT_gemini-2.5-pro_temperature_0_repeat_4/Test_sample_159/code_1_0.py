import math

def solve_spy_notebook_problem():
    """
    Calculates the highest number of observations James can record in his notebook.
    """
    # Step 1: Determine the number of states per page.
    # There are 3 basic fold types: Upper-right corner (U), Lower-right corner (L), and Vertical-half (V).
    # A page can have 0, 1, or at most 2 folds.

    # For the 80 middle pages, the folding order is not remembered.
    # - 0 folds: 1 state (unchanged)
    # - 1 fold: 3 states ({U}, {L}, {V})
    # - 2 folds: Combinations of 2 from 3 = 3 states ({U,L}, {U,V}, {L,V})
    states_order_doesnt_matter = 1 + 3 + 3
    num_pages_order_doesnt_matter = 80

    # For the first 10 and last 10 pages (20 total), the folding order is remembered.
    # - 0 folds: 1 state (unchanged)
    # - 1 fold: 3 states (U, L, V)
    # - 2 folds: The order can create distinct states.
    #   - Folding U then L is the same as L then U -> 1 state.
    #   - Folding U then V is different from V then U -> 2 states.
    #   - Folding L then V is different from V then L -> 2 states.
    #   - Total 2-fold states = 1 + 2 + 2 = 5 states.
    states_order_matters = 1 + 3 + 5
    num_pages_order_matters = 20

    # Step 2: Determine the number of unique observation types.
    # - 5 size estimates: few, small, medium, large, huge
    # - 6 time slots: 12am, 4am, 8am, 12pm, 4pm, 8pm
    num_sizes = 5
    num_times = 6
    observation_types = num_sizes * num_times

    # Step 3: Formulate the inequality and solve for N using logarithms.
    # The total number of notebook configurations must be >= the number of possible observation sequences.
    # (states_order_matters^20) * (states_order_doesnt_matter^80) >= observation_types^N
    # N <= (20 * log(9) + 80 * log(7)) / log(30)

    # Perform the calculation
    log_states_order_matters = math.log(states_order_matters)
    log_states_order_doesnt_matter = math.log(states_order_doesnt_matter)
    log_observation_types = math.log(observation_types)

    numerator = num_pages_order_matters * log_states_order_matters + num_pages_order_doesnt_matter * log_states_order_doesnt_matter
    denominator = log_observation_types

    max_n_float = numerator / denominator
    max_n_integer = math.floor(max_n_float)

    # Step 4: Print the explanation and the result.
    print("To find the highest number of observations (N), we solve the inequality:")
    print(f"({states_order_matters}^{num_pages_order_matters}) * ({states_order_doesnt_matter}^{num_pages_order_doesnt_matter}) >= {observation_types}^N")
    print("\nUsing logarithms, the equation to find the maximum N is:")
    print(f"N <= ({num_pages_order_matters} * log({states_order_matters}) + {num_pages_order_doesnt_matter} * log({states_order_doesnt_matter})) / log({observation_types})")
    
    print("\nPlugging in the numbers for the final equation:")
    final_equation_str = (
        f"floor( ("
        f"{num_pages_order_matters} * {log_states_order_matters:.4f} + "
        f"{num_pages_order_doesnt_matter} * {log_states_order_doesnt_matter:.4f}"
        f") / {log_observation_types:.4f} )"
    )
    print(final_equation_str)
    
    final_calculation_str = f"= floor( ({numerator:.4f}) / {denominator:.4f} ) = floor({max_n_float:.4f}) = {max_n_integer}"
    print(final_calculation_str)

    print(f"\nThe highest number of observations James can record is {max_n_integer}.")

solve_spy_notebook_problem()
<<<58>>>