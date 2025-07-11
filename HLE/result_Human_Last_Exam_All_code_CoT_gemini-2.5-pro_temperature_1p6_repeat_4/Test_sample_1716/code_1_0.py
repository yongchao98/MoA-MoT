import math

def solve_utility():
    """
    Calculates Alice's expected utility based on the principle of superrationality.
    """
    
    # Step 1 & 2: Define payoffs for symmetric outcomes, based on the problem and the principle of superrationality.
    # The problem states Alice and Bob are "superrational", which means they know they will reason
    # identically and thus make the same choice. We only need to consider the payoffs where
    # both choose the same action.
    
    # If both rest, their payoff is 0 each.
    payoff_rest = 0
    
    # If both bike, their payoff is -2 each.
    payoff_bike = -2
    
    # If both run, their payoff is -3 each.
    payoff_run = -3
    
    # Step 3 & 4: As rational agents, they will analyze these symmetric outcomes and choose
    # the one that gives the highest payoff.
    symmetric_payoffs = [payoff_rest, payoff_bike, payoff_run]
    expected_utility = max(symmetric_payoffs)

    print("This is a game of strategy between two 'superrational' agents, Alice and Bob.")
    print("Superrationality implies that both agents will follow the same logic and arrive at the same choice.")
    print("Therefore, we only need to consider the scenarios where both choose the same action:\n")
    print(f"1. If both Alice and Bob Rest, Alice's payoff is: {payoff_rest}")
    print(f"2. If both Alice and Bob Bike, Alice's payoff is: {payoff_bike}")
    print(f"3. If both Alice and Bob Run, Alice's payoff is: {payoff_run}\n")

    print("To maximize her own utility, Alice will choose the action from this set that results in the highest payoff.")
    print("Bob will reason in the exact same way.\n")
    
    # Final step: Print the calculation for Alice's expected utility.
    print("Alice's expected utility is the maximum of these possible outcomes.")
    print(f"Final Calculation: max({payoff_rest}, {payoff_bike}, {payoff_run}) = {expected_utility}")

solve_utility()