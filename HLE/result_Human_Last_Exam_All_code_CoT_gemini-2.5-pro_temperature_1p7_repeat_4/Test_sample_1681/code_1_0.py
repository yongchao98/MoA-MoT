def solve_quantum_dilemma():
    """
    Calculates the equilibrium payoff for the quantum Prisoner's Dilemma.
    """
    # Payoff values from the matrix for (C,C), (C,D), (D,C), and (D,D) outcomes
    payoff_CC = 5  # Reward for mutual cooperation
    payoff_CD_A = 0  # Sucker's payoff for player A
    payoff_DC_A = 7  # Temptation payoff for player A
    payoff_DD = 1  # Punishment for mutual defection

    # In the quantum game's (Q, Q) equilibrium, the probability of the
    # (Cooperate, Cooperate) outcome is 1, and all others are 0.
    prob_CC = 1.0
    prob_CD = 0.0
    prob_DC = 0.0
    prob_DD = 0.0

    # Calculate the expected payoff for a player
    equilibrium_payoff = (prob_CC * payoff_CC) + (prob_CD * payoff_CD_A) + (prob_DC * payoff_DC_A) + (prob_DD * payoff_DD)

    print("The quantum equilibrium is achieved when both players adopt the quantum strategy 'Q'.")
    print("This strategy results in the 'Cooperate, Cooperate' outcome with certainty.")
    print("\nThe payoff for each player at this equilibrium point is calculated based on the probabilities of each outcome:")
    
    # Print the equation with the numbers plugged in
    print(f"Payoff = P(CC) * Payoff(CC) + P(CD) * Payoff(CD) + P(DC) * Payoff(DC) + P(DD) * Payoff(DD)")
    print(f"Payoff = {prob_CC} * {payoff_CC} + {prob_CD} * {payoff_CD_A} + {prob_DC} * {payoff_DC_A} + {prob_DD} * {payoff_DD} = {int(equilibrium_payoff)}")


solve_quantum_dilemma()