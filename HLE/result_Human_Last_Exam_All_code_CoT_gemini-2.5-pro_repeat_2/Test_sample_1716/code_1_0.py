def calculate_expected_utility():
    """
    Calculates Alice's expected utility based on superrationality.

    In this game, Alice and Bob are superrational agents, meaning they are
    identical reasoners and know that they will both independently choose the
    same action. Therefore, we only need to consider the payoffs on the
    diagonal of the game matrix, where their actions are identical.

    The payoffs for Alice in these symmetric scenarios are:
    - If both Rest: 0
    - If both Bike: -2
    - If both Run: -3

    Alice, being rational, will choose the action that maximizes her payoff
    from these options.
    """
    
    # Payoffs for Alice when both she and Bob choose the same action
    payoff_rest_rest = 0
    payoff_bike_bike = -2
    payoff_run_run = -3
    
    # A superrational agent will choose the action that leads to the best symmetric outcome.
    possible_outcomes = [payoff_rest_rest, payoff_bike_bike, payoff_run_run]
    
    # Calculate Alice's utility, which is the maximum of the possible outcomes.
    alices_utility = max(possible_outcomes)
    
    # Print the equation as requested
    print(f"Alice's expected utility = max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {alices_utility}")

calculate_expected_utility()