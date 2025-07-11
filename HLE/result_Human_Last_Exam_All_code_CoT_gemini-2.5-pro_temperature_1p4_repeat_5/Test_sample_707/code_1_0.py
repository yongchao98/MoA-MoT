def agent_a_miscalculation():
    """
    This function demonstrates Agent A's flawed calculation.
    """
    initial_probability = 0.01
    quarterly_increase = 0.01
    total_quarters = 20

    # The goal is set at the start of the first quarter.
    # There are 19 subsequent quarters where an increase is added.
    num_increases = total_quarters - 1

    # Agent A's flawed calculation: linearly adding the increases
    # to the initial probability.
    final_probability = initial_probability + (num_increases * quarterly_increase)

    print("Agent A's flawed calculation for the cumulative probability of failure:")
    print(f"{final_probability:.2f} = {initial_probability} + {num_increases} * {quarterly_increase}")
    print("\nThis is incorrect because cumulative probability cannot be calculated by linearly adding the probability increases of each period.")

agent_a_miscalculation()