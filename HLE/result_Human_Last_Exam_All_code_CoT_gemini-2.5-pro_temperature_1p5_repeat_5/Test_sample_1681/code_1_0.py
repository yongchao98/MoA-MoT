import numpy as np

def calculate_quantum_equilibrium():
    """
    Calculates the payoffs for the Prisoner's Dilemma at the quantum Nash Equilibrium.

    In the quantum version of the game using the EWL protocol with maximal entanglement,
    a new equilibrium strategy Q emerges. When both players adopt Q, the outcome
    is mutual cooperation with probability 1.
    """
    # Payoff matrix values
    # R: Reward for mutual cooperation
    # S: Sucker's payoff for cooperating while the other defects
    # T: Temptation to defect while the other cooperates
    # P: Punishment for mutual defection
    R, S, T, P = 5, 0, 7, 1

    # At the quantum Nash Equilibrium (Q, Q), the outcome probabilities are fixed:
    P_CC = 1.0  # Probability of (Cooperate, Cooperate)
    P_CD = 0.0  # Probability of (Cooperate, Defect)
    P_DC = 0.0  # Probability of (Defect, Cooperate)
    P_DD = 0.0  # Probability of (Defect, Defect)

    # Calculate Alice's payoff (Player 1)
    # Payoff_A = R * P_CC + S * P_CD + T * P_DC + P * P_DD
    alice_payoff = R * P_CC + S * P_CD + T * P_DC + P * P_DD

    # Calculate Bob's payoff (Player 2)
    # The payoff for Bob depends on Alice's choice.
    # If Alice Cooperates (first outcome), Bob gets R or T.
    # If Alice Defects (first outcome), Bob gets S or P.
    # Payoff_B = R * P_CC + T * P_CD + S * P_DC + P * P_DD
    bob_payoff = R * P_CC + T * P_CD + S * P_DC + P * P_DD

    print("Quantum Prisoner's Dilemma Equilibrium Calculation")
    print("-" * 50)
    print("The equilibrium is achieved when both players adopt the quantum strategy Q.")
    print("This results in the following outcome probabilities:")
    print(f"P(Cooperate, Cooperate) = {P_CC}")
    print(f"P(Cooperate, Defect)   = {P_CD}")
    print(f"P(Defect, Cooperate)   = {P_DC}")
    print(f"P(Defect, Defect)     = {P_DD}")
    print("-" * 50)

    # Print the detailed payoff calculation for Alice
    print("Alice's payoff calculation:")
    print(f"Payoff_A = {R} * {P_CC} + {S} * {P_CD} + {T} * {P_DC} + {P} * {P_DD}")
    print(f"         = {alice_payoff}")

    # Print the detailed payoff calculation for Bob
    print("\nBob's payoff calculation:")
    print(f"Payoff_B = {R} * {P_CC} + {T} * {P_CD} + {S} * {P_DC} + {P} * {P_DD}")
    print(f"         = {bob_payoff}")
    
    print("\nThus, the equilibrium point yields a payoff of ({}, {})".format(int(alice_payoff), int(bob_payoff)))

calculate_quantum_equilibrium()
