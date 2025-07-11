import collections

def solve_alice_utility():
    """
    Solves for Alice's expected utility by modeling the game and applying
    the concept of superrationality to find the logical outcome.
    """
    # Define players' actions
    actions = ['rest', 'bike', 'run']

    # Define Alice's payoff matrix: U(Alice's Action, Bob's Action)
    payoff_alice = {
        'rest': {'rest': 0, 'bike': 2, 'run': 4},
        'bike': {'rest': 0, 'bike': -2, 'run': 2},
        'run':  {'rest': 0, 'bike': 0, 'run': -3}
    }

    print("This script calculates Alice's expected utility.")
    print("The problem describes a symmetric game between two 'superrational' agents, Alice and Bob.")
    print("A superrational agent will eliminate weakly dominated strategies to make a decision.\n")

    print("Step 1: Analyze Alice's strategies for dominance.")
    print("A strategy is weakly dominated if another strategy is always at least as good, and strictly better in at least one case.")
    print("-" * 50)

    undominated_strategies = list(actions)
    dominated_by = {}

    # Iterate through all pairs of strategies to check for dominance
    for s1 in actions:
        for s2 in actions:
            if s1 == s2:
                continue

            # Check if s1 weakly dominates s2
            is_at_least_as_good = all(payoff_alice[s1][opp_action] >= payoff_alice[s2][opp_action] for opp_action in actions)
            is_strictly_better = any(payoff_alice[s1][opp_action] > payoff_alice[s2][opp_action] for opp_action in actions)

            if is_at_least_as_good and is_strictly_better:
                # s1 dominates s2
                print(f"Analysis: Strategy '{s1}' weakly dominates strategy '{s2}'.")
                print(f"  - Payoffs for '{s1}': (vs Rest: {payoff_alice[s1]['rest']}, vs Bike: {payoff_alice[s1]['bike']}, vs Run: {payoff_alice[s1]['run']})")
                print(f"  - Payoffs for '{s2}': (vs Rest: {payoff_alice[s2]['rest']}, vs Bike: {payoff_alice[s2]['bike']}, vs Run: {payof_alice[s2]['run']})")
                print("")
                if s2 in undominated_strategies:
                    undominated_strategies.remove(s2)
                dominated_by[s2] = s1

    print("-" * 50)
    print("Step 2: Determine the rational choice for Alice and Bob.")

    if len(undominated_strategies) == 1:
        alice_choice = undominated_strategies[0]
        print(f"After eliminating the dominated strategies ('bike' and 'run'), the only undominated strategy for Alice is '{alice_choice}'.")
        # Due to symmetry and superrationality, Bob makes the same choice.
        bob_choice = alice_choice
        print(f"As a superrational agent, Bob will reason identically and also choose '{bob_choice}'.\n")
    else:
        print("Could not determine a single rational choice through dominance.")
        return

    print("Step 3: Calculate Alice's expected utility from the outcome.")
    final_payoff = payoff_alice[alice_choice][bob_choice]
    print(f"The determined outcome of the game is (Alice: {alice_choice}, Bob: {bob_choice}).")
    
    print("\nThe final equation for Alice's expected utility is based on this single, determined outcome:")
    print(f"Alice's Expected Utility = Payoff(Alice='{alice_choice}', Bob='{bob_choice}')")
    print(f"Alice's Expected Utility = {final_payoff}")

solve_alice_utility()