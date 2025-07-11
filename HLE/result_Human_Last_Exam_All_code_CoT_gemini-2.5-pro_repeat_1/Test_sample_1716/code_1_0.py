import numpy as np

def solve_utility():
    """
    Calculates Alice's expected utility based on superrationality and spite.
    """

    # Step 1: Model the game with a payoff matrix (Alice's payoff, Bob's payoff)
    # Rows: Alice's actions (Rest, Bike, Run)
    # Columns: Bob's actions (Rest, Bike, Run)
    payoff_matrix = {
        'Rest': {'Rest': (0, 0), 'Bike': (2, 0), 'Run': (4, 0)},
        'Bike': {'Rest': (0, 2), 'Bike': (-2, -2), 'Run': (2, 0)},
        'Run':  {'Rest': (0, 4), 'Bike': (0, 2), 'Run': (-3, -3)}
    }
    actions = ['Rest', 'Bike', 'Run']
    alice_payoffs = np.array([[payoff_matrix[r][c][0] for c in actions] for r in actions])
    bob_payoffs = np.array([[payoff_matrix[r][c][1] for c in actions] for r in actions])

    print("Step 1: The payoff matrix (Alice's Payoff, Bob's Payoff) is:")
    for i, row_action in enumerate(actions):
        row_str = f"Alice {row_action}:".ljust(12)
        for j, col_action in enumerate(actions):
            row_str += f" vs Bob {col_action}: ({alice_payoffs[i,j]}, {bob_payoffs[i,j]})".ljust(25)
        print(row_str)
    print("-" * 80)

    # Step 2 & 3: Alice's superrational decision
    # Alice's 'Rest' strategy payoffs: (0, 2, 4)
    # Alice's 'Bike' strategy payoffs: (0, -2, 2)
    # Alice's 'Run'  strategy payoffs: (0, 0, -3)
    # Rest weakly dominates Bike because 0>=0, 2>-2, 4>2.
    # Rest weakly dominates Run because 0>=0, 2>0, 4>-3.
    print("Step 2: Analyzing Alice's strategies.")
    print("Alice's payoffs if she Rests are (0, 2, 4) depending on Bob's action.")
    print("Alice's payoffs if she Bikes are (0, -2, 2).")
    print("Alice's payoffs if she Runs are (0, 0, -3).")
    print("\nThe 'Rest' strategy is always at least as good as, and sometimes better than, her other options.")
    print("This means 'Rest' weakly dominates 'Bike' and 'Run'.")
    print("\nStep 3: As a superrational agent, Alice will eliminate her weakly dominated strategies and choose to 'Rest'.")
    alice_choice = 'Rest'
    print(f"Alice's determined action: {alice_choice}")
    print("-" * 80)

    # Step 4 & 5: Bob's superrational prediction and decision
    print("Step 4: Bob, also being superrational, simulates Alice's reasoning and predicts she will 'Rest'.")
    # Get Bob's payoffs when Alice rests (the first row of Bob's payoff matrix)
    bob_payoffs_vs_alice_rest = bob_payoffs[0, :]
    print(f"Bob analyzes his payoffs assuming Alice Rests: for his actions {{Rest, Bike, Run}}, his payoffs are {tuple(bob_payoffs_vs_alice_rest)}.")
    print("Bob sees his payoff is 0 no matter what he does. He is indifferent.")
    print("\nStep 5: To break the tie, Bob uses the fact that he despises Alice.")
    print("He will choose the action that minimizes Alice's payoff.")

    # Get Alice's payoffs when she rests
    alice_payoffs_if_resting = alice_payoffs[0, :]
    print(f"Alice's payoffs when she rests are {tuple(alice_payoffs_if_resting)} for Bob's actions {{Rest, Bike, Run}}.")

    # Find the action Bob takes to minimize Alice's payoff
    min_alice_payoff = min(alice_payoffs_if_resting)
    bob_choice_index = list(alice_payoffs_if_resting).index(min_alice_payoff)
    bob_choice = actions[bob_choice_index]
    print(f"The minimum payoff for Alice is {min_alice_payoff}, which occurs when Bob chooses '{bob_choice}'.")
    print(f"Bob's determined action: {bob_choice}")
    print("-" * 80)

    # Step 6 & 7: Final Outcome and Alice's Utility
    print("Step 6: The determined outcome is (Alice: Rest, Bob: Rest).")
    final_payoff = payoff_matrix[alice_choice][bob_choice]
    alices_final_utility = final_payoff[0]

    print("\nStep 7: Alice's expected utility is her payoff in this outcome.")
    print(f"Alice's Expected Utility = Payoff({alice_choice}, {bob_choice}) = {alices_final_utility}")
    # The final equation is trivial, but we print it as requested.
    print(f"\nFinal Equation:")
    print(f"Alice's Utility = {alices_final_utility}")

    return alices_final_utility

if __name__ == '__main__':
    utility = solve_utility()
    # The final answer format as requested by the user.
    # print(f"\n<<<{utility}>>>")