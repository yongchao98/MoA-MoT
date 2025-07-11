import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma.
    """
    # Payoff Matrix values: R (Reward), S (Sucker), T (Temptation), P (Punishment)
    # Corresponding to (C,C), (C,D), (D,C), (D,D)
    # Player A is the row player, Player B is the column player.
    # Payoff for A: (R, S, T, P) = (5, 0, 7, 1)
    # Payoff for B: (R, T, S, P) = (5, 7, 0, 1)
    R, S, T, P = 5, 0, 7, 1

    # Define standard 2x2 matrices (quantum gates)
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- Define Player Strategies as Unitary Matrices ---
    # Cooperate (C) is the identity operator
    C_gate = I
    # Defect (D) is a bit-flip. The 'i' makes it unitary with determinant 1.
    D_gate = 1j * sigma_x
    # The special Quantum strategy (Q) that resolves the dilemma
    Q_gate = 1j * sigma_z

    # --- Setup the Quantum Game (EWL Protocol) ---
    # Initial state is |00>
    psi_initial = np.array([1, 0, 0, 0], dtype=complex)

    # Entangling operator J
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    # Disentangling operator J_dag (conjugate transpose)
    J_dag = J.conj().T

    # --- Analyze the Equilibrium Point (Q, Q) ---
    # Both players choose the quantum strategy Q
    U_A = Q_gate
    U_B = Q_gate

    # The combined strategy operator for the two players
    U_total = np.kron(U_A, U_B)

    # Calculate the final state before measurement
    psi_final = J_dag @ U_total @ J @ psi_initial

    # Probabilities are the squared magnitudes of the final state vector's components
    # The components correspond to |00>, |01>, |10>, |11>
    probabilities = np.abs(psi_final)**2
    p_cc, p_cd, p_dc, p_dd = probabilities[0], probabilities[1], probabilities[2], probabilities[3]

    # Calculate the expected payoff for each player
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd

    # --- Print the Result ---
    # When both players choose the quantum strategy Q, a new, stable equilibrium is reached.
    # The resulting outcome corresponds to the (Cooperate, Cooperate) case.
    print("Quantum Equilibrium Point (Q, Q) Analysis:")
    print(f"Probabilities of classical outcomes: P(CC)={p_cc:.2f}, P(CD)={p_cd:.2f}, P(DC)={p_dc:.2f}, P(DD)={p_dd:.2f}")
    print("\nExpected Payoff for Player A:")
    print(f"Equation: ({R} * {p_cc:.2f}) + ({S} * {p_cd:.2f}) + ({T} * {p_dc:.2f}) + ({P} * {p_dd:.2f}) = {payoff_A:.2f}")
    
    print("\nExpected Payoff for Player B:")
    print(f"Equation: ({R} * {p_cc:.2f}) + ({T} * {p_cd:.2f}) + ({S} * {p_dc:.2f}) + ({P} * {p_dd:.2f}) = {payoff_B:.2f}")
    
    print(f"\nThe quantum equilibrium point results in the payoff vector ({payoff_A:.1f}, {payoff_B:.1f}).")

solve_quantum_prisoners_dilemma()
<<<5>>>