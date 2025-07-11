def solve_utility():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    # Step 1: Define the payoffs for the symmetric outcomes based on the problem description.
    # If both rest, payoff is 0.
    # If both bike, payoff is -2.
    # If both run, payoff is -3.
    payoff_rest_rest = 0
    payoff_bike_bike = -2
    payoff_run_run = -3

    # Step 2: Explain the reasoning based on superrationality.
    print("Alice and Bob are superrational agents in a symmetric game.")
    print("This means they will both independently deduce and select the same optimal strategy.")
    print("Therefore, we only need to consider the outcomes where they both choose the same action:\n")

    # Step 3: List the payoffs for the relevant outcomes.
    print(f"1. If both Alice and Bob Rest, Alice's payoff is {payoff_rest_rest}.")
    print(f"2. If both Alice and Bob Bike, Alice's payoff is {payoff_bike_bike}.")
    print(f"3. If both Alice and Bob Run, Alice's payoff is {payoff_run_run}.\n")

    # Step 4: Determine the optimal choice by finding the maximum payoff.
    # A rational agent will choose the action that leads to the best outcome for themselves.
    alices_expected_utility = max(payoff_rest_rest, payoff_bike_bike, payoff_run_run)

    # Step 5: Print the final calculation and conclusion.
    print("To maximize her utility, Alice will choose the action that results in the highest payoff among these symmetric outcomes.")
    print(f"The final equation is: max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run})")
    print(f"Alice's expected utility is the result of this equation, which is {alices_expected_utility}.")

solve_utility()
<<<0>>>