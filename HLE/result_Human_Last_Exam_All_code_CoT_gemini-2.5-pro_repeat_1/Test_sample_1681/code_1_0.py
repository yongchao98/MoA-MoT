import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma.

    In the EWL quantum game framework with maximal entanglement, a new Nash
    Equilibrium emerges where all four classical outcomes become equally likely.
    The equilibrium payoff is the average of all possible payoffs for a player.
    """

    # Payoff matrix:
    # rows are Alice's choices (C, D), columns are Bob's choices (C, D)
    # payoffs are (Alice's payoff, Bob's payoff)
    payoffs = {
        'CC': (5, 5),
        'CD': (0, 7),
        'DC': (7, 0),
        'DD': (1, 1)
    }

    # At the quantum Nash Equilibrium, the probability of each outcome is 0.25
    prob_CC = 0.25
    prob_CD = 0.25
    prob_DC = 0.25
    prob_DD = 0.25

    # Alice's payoffs for each outcome
    alice_payoff_CC = payoffs['CC'][0]
    alice_payoff_CD = payoffs['CD'][0]
    alice_payoff_DC = payoffs['DC'][0]
    alice_payoff_DD = payoffs['DD'][0]

    # Calculate Alice's expected payoff at the equilibrium
    alice_expected_payoff = (
        alice_payoff_CC * prob_CC +
        alice_payoff_CD * prob_CD +
        alice_payoff_DC * prob_DC +
        alice_payoff_DD * prob_DD
    )

    # Since the game is symmetric, Bob's expected payoff is the same.
    # We can calculate it for confirmation.
    bob_payoff_CC = payoffs['CC'][1]
    bob_payoff_CD = payoffs['CD'][1]
    bob_payoff_DC = payoffs['DC'][1]
    bob_payoff_DD = payoffs['DD'][1]

    bob_expected_payoff = (
        bob_payoff_CC * prob_CC +
        bob_payoff_CD * prob_CD +
        bob_payoff_DC * prob_DC +
        bob_payoff_DD * prob_DD
    )
    
    # The equilibrium point is the payoff pair (Alice's payoff, Bob's payoff)
    # The question asks for "the" equilibrium point, which is the value of the symmetric payoff.
    equilibrium_point = alice_expected_payoff

    print("Quantum Prisoner's Dilemma Equilibrium Calculation:")
    print(f"The payoff matrix provides the following values for Alice:")
    print(f"  (Cooperate, Cooperate): {alice_payoff_CC}")
    print(f"  (Cooperate, Defect):   {alice_payoff_CD}")
    print(f"  (Defect,   Cooperate): {alice_payoff_DC}")
    print(f"  (Defect,   Defect):   {alice_payoff_DD}")
    print("\nAt the quantum equilibrium, each outcome has a probability of 0.25.")
    print("The expected payoff is the sum of (payoff * probability) for all outcomes.")
    print("\nFinal Payoff Equation:")
    print(f"Payoff = ({alice_payoff_CC} * {prob_CC}) + ({alice_payoff_CD} * {prob_CD}) + ({alice_payoff_DC} * {prob_DC}) + ({alice_payoff_DD} * {prob_DD})")
    print(f"Payoff = {alice_payoff_CC * prob_CC} + {alice_payoff_CD * prob_CD} + {alice_payoff_DC * prob_DC} + {alice_payoff_DD * prob_DD}")
    print(f"Equilibrium Payoff = {equilibrium_point}")
    print(f"\nThe equilibrium point is a symmetric payoff of ({alice_expected_payoff}, {bob_expected_payoff}).")

solve_quantum_prisoners_dilemma()
<<<3.25>>>