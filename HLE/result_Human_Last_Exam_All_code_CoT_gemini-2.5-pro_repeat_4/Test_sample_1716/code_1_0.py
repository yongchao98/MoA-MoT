def solve_utility():
    """
    Calculates Alice's expected utility in a game with a superrational opponent.
    """
    # Payoffs for Alice are defined based on the problem description.
    # (Alice's Action, Bob's Action) -> Alice's Payoff
    # This can be represented as a matrix where rows are Alice's choices
    # and columns are Bob's choices.

    # Payoff when Alice Rests and Bob Runs:
    payoff_alice_rests_bob_runs = 4

    # Payoff when Alice Runs and Bob Rests:
    payoff_alice_runs_bob_rests = 0
    
    # The problem states Alice and Bob are superrational. In a symmetric game,
    # they will coordinate on the symmetric outcome that maximizes their mutual payoff.
    # The pure symmetric outcomes are (Rest,Rest)->0, (Bike,Bike)->-2, (Run,Run)->-3.
    # A better symmetric outcome can be achieved by coordinating on a 50/50 probability
    # of playing (Alice Rests, Bob Runs) and (Alice Runs, Bob Rests).
    # This yields the highest possible symmetric payoff for both players.
    
    # Probabilities for the coordinated strategy
    prob_one = 0.5
    prob_two = 0.5

    # Calculate Alice's expected utility from this optimal strategy
    expected_utility = (prob_one * payoff_alice_rests_bob_runs) + (prob_two * payoff_alice_runs_bob_rests)

    # Print the final equation showing the calculation
    print(f"Alice's expected utility is derived from the optimal symmetric strategy:")
    print(f"({prob_one} * {payoff_alice_rests_bob_runs}) + ({prob_two} * {payoff_alice_runs_bob_rests}) = {expected_utility}")

solve_utility()