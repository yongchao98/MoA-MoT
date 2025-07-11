import sys

def solve_alice_utility():
    """
    Calculates Alice's expected utility based on superrationality in a symmetric game.
    """
    # Step 1: Define the payoff matrix for Alice based on the problem description.
    # The matrix represents Payoff(Alice's Action, Bob's Action).
    # Rows: Alice's choice (Rest, Bike, Run)
    # Columns: Bob's choice (Rest, Bike, Run)
    payoffs = {
        "Rest": {"Rest": 0, "Bike": 2, "Run": 4},
        "Bike": {"Rest": 0, "Bike": -2, "Run": 2},
        "Run":  {"Rest": 0, "Bike": 0, "Run": -3},
    }

    print("This is a symmetric game between two superrational agents, Alice and Bob.")
    print("The payoff matrix for Alice is as follows:")
    print("       Bob's Action")
    print("Alice's   | Rest | Bike | Run")
    print("----------------------------")
    print(f"Rest      | {payoffs['Rest']['Rest']:^4} | {payoffs['Rest']['Bike']:^4} | {payoffs['Rest']['Run']:^4}")
    print(f"Bike      | {payoffs['Bike']['Rest']:^4} | {payoffs['Bike']['Bike']:^4} | {payoffs['Bike']['Run']:^4}")
    print(f"Run       | {payoffs['Run']['Rest']:^4} | {payoffs['Run']['Bike']:^4} | {payoffs['Run']['Run']:^4}")
    print("\n")

    # Step 2: Apply the principle of superrationality.
    # Since Alice and Bob are identically rational and self-aware, they will independently
    # deduce and commit to the same strategy. This means we only need to evaluate
    # the symmetric outcomes where they both perform the same action.
    print("Reasoning:")
    print("1. As superrational agents in a symmetric game, both Alice and Bob know they will inevitably arrive at the same conclusion and choose the same action.")
    print("2. Therefore, Alice must evaluate her choices assuming Bob will make the same choice she does.")
    print("3. Alice's potential payoffs are limited to the diagonal of the matrix:")
    
    payoff_rest_rest = payoffs["Rest"]["Rest"]
    payoff_bike_bike = payoffs["Bike"]["Bike"]
    payoff_run_run = payoffs["Run"]["Run"]
    
    print(f"   - If both Rest: Alice's payoff is {payoff_rest_rest}")
    print(f"   - If both Bike: Alice's payoff is {payoff_bike_bike}")
    print(f"   - If both Run:  Alice's payoff is {payoff_run_run}")
    print("\n")

    # Step 3: Determine the optimal choice.
    # A rational agent maximizes their own utility. Alice will choose the action that
    # yields the highest payoff among the symmetric outcomes.
    print("To maximize her own utility, Alice compares these payoffs and chooses the action associated with the highest value.")
    print(f"Comparing {payoff_rest_rest}, {payoff_bike_bike}, and {payoff_run_run}, the best outcome for her is {max(payoff_rest_rest, payoff_bike_bike, payoff_run_run)}.")
    print("This corresponds to both agents choosing to Rest.")
    print("\n")

    # Step 4: Calculate and display the final expected utility.
    # Because the choice is deterministic, Alice's expected utility is simply her payoff
    # in the single resulting scenario: (Rest, Rest).
    # We can write this as an equation where the probability of choosing Rest is 1, and 0 for other actions.
    
    # Probabilities of choosing the symmetric outcome (Rest, Bike, or Run)
    prob_rest = 1
    prob_bike = 0
    prob_run = 0

    expected_utility = (prob_rest * payoff_rest_rest) + \
                       (prob_bike * payoff_bike_bike) + \
                       (prob_run * payoff_run_run)
                       
    print("Alice's expected utility is calculated based on the probabilities of each symmetric outcome occurring.")
    print("Since 'Rest' is the optimal choice, its probability is 1, while the others are 0.")
    print("Final Calculation:")
    print(f"E[Utility] = P(Rest) * Payoff(Rest,Rest) + P(Bike) * Payoff(Bike,Bike) + P(Run) * Payoff(Run,Run)")
    print(f"E[Utility] = {prob_rest} * {payoff_rest_rest} + {prob_bike} * {payoff_bike_bike} + {prob_run} * {payoff_run_run} = {expected_utility}")
    
    # Appending final answer in the required format to the stream
    sys.stdout.write("\n<<<0>>>")

solve_alice_utility()