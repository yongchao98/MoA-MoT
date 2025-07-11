def solve_quantum_dilemma():
    """
    Calculates and analyzes the equilibrium point for the quantum prisoner's dilemma.
    """
    # Step 1: Define the classical payoff matrix values
    # (Player1_payoff, Player2_payoff)
    # R: (Cooperate, Cooperate)
    # S: (Cooperate, Defect) for player 1
    # T: (Defect, Cooperate) for player 1
    # P: (Defect, Defect)
    R = 5
    S = 0
    T = 7
    P = 1

    print("Original Payoff Matrix (Player 1, Player 2):")
    print(f"           Player 2 Cooperate    Player 2 Defect")
    print(f"Player 1 Cooperate: ({R}, {R})             ({S}, {T})")
    print(f"Player 1 Defect:   ({T}, {S})             ({P}, {P})")
    print("-" * 50)

    # Step 2: Under the EWL quantum game protocol with maximal entanglement,
    # the game is equivalent to a classical game with transformed payoffs.
    print("Calculating payoffs for the equivalent game under the EWL quantum protocol...")

    R_q = R
    P_q = P
    S_q = (S + T) / 2
    T_q = (T + S) / 2 # Note: S_q will equal T_q

    print(f"The 'Cooperate, Cooperate' payoff R' remains R: {R_q}")
    print(f"The 'Defect, Defect' payoff P' remains P: {P_q}")
    print(f"The 'Cooperate, Defect' payoff S' becomes (S + T) / 2: ({S} + {T}) / 2 = {S_q}")
    print(f"The 'Defect, Cooperate' payoff T' becomes (T + S) / 2: ({T} + {S}) / 2 = {T_q}")
    print("-" * 50)

    # Step 3: Construct and print the new game matrix
    print("New Equivalent Payoff Matrix (Player 1, Player 2):")
    print(f"           Player 2 Cooperate    Player 2 Defect")
    print(f"Player 1 Cooperate: ({R_q}, {R_q})           ({S_q}, {T_q})")
    print(f"Player 1 Defect:   ({T_q}, {S_q})           ({P_q}, {P_q})")
    print("-" * 50)


    # Step 4: Analyze the new game for dominant strategies
    print("Analyzing Player 1's strategy:")
    # Case 1: Player 2 Cooperates
    print("If Player 2 Cooperates, Player 1's payoffs are:")
    print(f"  - Cooperate: {R_q}")
    print(f"  - Defect: {T_q}")
    print(f"Player 1 should Cooperate because {R_q} > {T_q}")
    print("-" * 25)

    # Case 2: Player 2 Defects
    print("If Player 2 Defects, Player 1's payoffs are:")
    print(f"  - Cooperate: {S_q}")
    print(f"  - Defect: {P_q}")
    print(f"Player 1 should Cooperate because {S_q} > {P_q}")
    print("-" * 25)

    print("\nConclusion: 'Cooperate' is the dominant strategy for Player 1.")
    print("Due to the symmetry of the game, 'Cooperate' is also the dominant strategy for Player 2.")
    print("-" * 50)

    # Step 5: State the equilibrium
    print(f"The unique Nash Equilibrium is (Cooperate, Cooperate).")
    print(f"The payoff at this equilibrium is ({R_q}, {R_q}).")
    print("\nThe dilemma is resolved: the mutually beneficial outcome is now a stable equilibrium.")


if __name__ == '__main__':
    solve_quantum_dilemma()