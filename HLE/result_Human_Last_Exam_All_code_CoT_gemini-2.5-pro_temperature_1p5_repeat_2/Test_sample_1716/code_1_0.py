def solve_utility():
    """
    Calculates Alice's expected utility based on the principles of superrationality.
    """
    # Step 1: Define the payoffs for the symmetric outcomes based on the problem description.
    # In a superrational game, agents know they will make the same choice.
    # We only need to evaluate the scenarios where Alice and Bob do the same thing.
    
    # Payoff for Alice if both rest
    payoff_rest_rest = 0
    # Payoff for Alice if both bike
    payoff_bike_bike = -2
    # Payoff for Alice if both run
    payoff_run_run = -3

    # Step 2: Explain the reasoning.
    print("Superrationality assumes that identical, perfectly rational agents will choose the same strategy.")
    print("Therefore, we only need to compare the outcomes where Alice and Bob both make the same choice:")
    print(f"- If both rest, Alice's utility is {payoff_rest_rest}.")
    print(f"- If both bike, Alice's utility is {payoff_bike_bike}.")
    print(f"- If both run, Alice's utility is {payoff_run_run}.")
    print("\n")
    
    # Step 3: Alice chooses the action that leads to the best outcome for herself.
    # This is equivalent to finding the maximum of the possible payoffs.
    alices_expected_utility = max(payoff_rest_rest, payoff_bike_bike, payoff_run_run)

    # Step 4: Print the final calculation and the result.
    # The user wants the final equation printed with the numbers.
    print("To find her expected utility, Alice will choose the action that results in the maximum payoff.")
    print("The calculation is:")
    print(f"max({payoff_rest_rest}, {payoff_bike_bike}, {payoff_run_run}) = {alices_expected_utility}")
    print("\n")
    print(f"Thus, Alice's expected utility is {alices_expected_utility}.")

solve_utility()