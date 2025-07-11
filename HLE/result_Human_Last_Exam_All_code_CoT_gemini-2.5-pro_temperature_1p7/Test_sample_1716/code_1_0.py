def solve_alice_utility():
    """
    Calculates Alice's expected utility based on superrationality.
    """
    # Payoffs for Alice based on the actions of (Alice, Bob)
    payoffs = {
        ('Rest', 'Rest'): 0,
        ('Bike', 'Bike'): -2,
        ('Run', 'Run'): -3,
    }

    print("Step 1: Understand Superrationality")
    print("Alice and Bob are superrational. This means they know that since they are both perfectly logical, they will independently reason their way to the same conclusion.")
    print("Therefore, they will both choose the same action. We only need to analyze the outcomes where their choices are identical.")
    print("\n------------------------------------\n")

    print("Step 2: Identify Payoffs for Symmetric Outcomes")
    print("We look at Alice's payoff for each of the three symmetric scenarios:")
    
    payoff_rest = payoffs[('Rest', 'Rest')]
    payoff_bike = payoffs[('Bike', 'Bike')]
    payoff_run = payoffs[('Run', 'Run')]

    print(f"If both Alice and Bob Rest, Alice's payoff is: {payoff_rest}")
    print(f"If both Alice and Bob Bike, Alice's payoff is: {payoff_bike}")
    print(f"If both Alice and Bob Run, Alice's payoff is: {payoff_run}")
    print("\n------------------------------------\n")

    print("Step 3: Determine the Optimal Choice")
    print("To maximize her own utility, Alice will choose the action that leads to the highest payoff among these possibilities.")
    print("The decision is based on the following calculation:")
    
    # Find the maximum payoff from the symmetric outcomes
    optimal_payoff = max(payoff_rest, payoff_bike, payoff_run)
    
    print(f"max({payoff_rest}, {payoff_bike}, {payoff_run})")
    print("\n------------------------------------\n")

    print("Step 4: Final Expected Utility")
    print("The superrational choice is 'Rest', as it yields the highest payoff of 0.")
    print("Since this outcome is certain under superrationality, Alice's expected utility is equal to this payoff.")
    print(f"\nAlice's expected utility is: {optimal_payoff}")

solve_alice_utility()