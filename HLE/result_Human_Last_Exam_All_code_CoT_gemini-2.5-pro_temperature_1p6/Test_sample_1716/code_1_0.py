import sys

def solve_alice_utility():
    """
    Calculates Alice's expected utility based on the game theory puzzle.
    """
    # Step 1: Define Alice's payoff matrix based on the problem description.
    # The matrix is represented as a dictionary where keys are Alice's actions
    # and values are lists of payoffs corresponding to Bob's actions (Rest, Bike, Run).
    # Format: payoffs_alice[Alice_Action] = [vs_Bob_Rest, vs_Bob_Bike, vs_Bob_Run]
    payoffs_alice = {
        'Rest': [0, 2, 4],
        'Bike': [0, -2, 2],
        'Run':  [0, 0, -3]
    }
    actions = ['Rest', 'Bike', 'Run']
    action_indices = {action: i for i, action in enumerate(actions)}

    # Step 2 & 3: Apply the logic of superrationality.
    # Superrational agents know they will reason identically and make the same choice.
    # Therefore, we only need to evaluate the diagonal of the payoff matrix
    # where both choose the same action.
    payoff_rr = payoffs_alice['Rest'][action_indices['Rest']]
    payoff_bb = payoffs_alice['Bike'][action_indices['Bike']]
    payoff_kk = payoffs_alice['Run'][action_indices['Run']]

    diagonal_payoffs = [payoff_rr, payoff_bb, payoff_kk]

    # Step 4 & 5: Determine the optimal choice and the resulting utility.
    # Alice will choose the action that maximizes her utility among the symmetric outcomes.
    # Her expected utility is the result of this deterministic choice.
    alice_utility = max(diagonal_payoffs)

    # Print the explanation and final calculation as requested.
    print("Based on the problem description, Alice's payoff matrix is as follows (Rows: Alice's actions, Columns: Bob's actions):")
    print("      Bob: Rest  Bike   Run")
    print(f"Alice Rest:   {payoffs_alice['Rest'][0]:>2d}     {payoffs_alice['Rest'][1]:>2d}      {payoffs_alice['Rest'][2]:>2d}")
    print(f"Alice Bike:   {payoffs_alice['Bike'][0]:>2d}    {payoffs_alice['Bike'][1]:>2d}      {payoffs_alice['Bike'][2]:>2d}")
    print(f"Alice Run:    {payoffs_alice['Run'][0]:>2d}     {payoffs_alice['Run'][1]:>2d}     {payoffs_alice['Run'][2]:>2d}")

    print("\nAlice and Bob are superrational, so they know they will independently arrive at the same strategy.")
    print("This means we only need to consider the outcomes where they both perform the same action:")
    print(f"- If both Rest, Alice's utility is {payoff_rr}.")
    print(f"- If both Bike, Alice's utility is {payoff_bb}.")
    print(f"- If both Run,  Alice's utility is {payoff_kk}.")
    
    print("\nTo maximize her utility, Alice will choose the action that results in the highest payoff from these options.")
    print("Her expected utility is the maximum of these values.")
    print("\nFinal Calculation:")
    print(f"Alice's Expected Utility = max({payoff_rr}, {payoff_bb}, {payoff_kk})")
    print(f"Alice's Expected Utility = {alice_utility}")


solve_alice_utility()