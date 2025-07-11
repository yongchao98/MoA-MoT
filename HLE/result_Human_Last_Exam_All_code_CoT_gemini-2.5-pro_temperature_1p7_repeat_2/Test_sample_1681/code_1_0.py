import numpy as np

def run_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma.
    """
    # 1. Define Payoff Matrix values
    # Payoffs are (Player A, Player B) for outcomes (CC, CD, DC, DD)
    # R: Reward for mutual cooperation
    # S: Sucker's payoff
    # T: Temptation to defect
    # P: Punishment for mutual defection
    R, T, S, P = 5, 7, 0, 1
    payoffs_A = np.array([R, S, T, P])

    # 2. Define Quantum states and Operators
    # Basis vectors |C> = |0> and |D> = |1>
    C = np.array([1, 0])
    D = np.array([0, 1])

    # Basis for the 2-qubit system
    CC = np.kron(C, C)
    CD = np.kron(C, D)
    DC = np.kron(D, C)
    DD = np.kron(D, D)

    # Pauli Matrices
    I = np.eye(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # 3. Define Game Operators
    # Entangling Operator for maximal entanglement (gamma = pi/2)
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dagger = J.conj().T

    # Strategy Operators
    C_op = I                                 # Cooperate
    D_op = 1j * sigma_x                      # Defect
    Q_op = 1j * sigma_z                      # Quantum strategy

    # 4. Set player strategies to the quantum equilibrium
    U_A = Q_op
    U_B = Q_op

    # 5. Simulate the game
    # Initial state
    psi_initial = CC

    # Apply entangling operator
    psi_entangled = J @ psi_initial

    # Players apply their strategies
    players_op = np.kron(U_A, U_B)
    psi_moved = players_op @ psi_entangled

    # Apply disentangling operator before measurement
    psi_final = J_dagger @ psi_moved

    # 6. Calculate Probabilities
    # Probabilities are the squared magnitudes of the final state's amplitudes
    # projected onto the classical basis states.
    P_CC = np.abs(CC.conj().T @ psi_final)**2
    P_CD = np.abs(CD.conj().T @ psi_final)**2
    P_DC = np.abs(DC.conj().T @ psi_final)**2
    P_DD = np.abs(DD.conj().T @ psi_final)**2

    # 7. Calculate Equilibrium Payoff
    player_A_payoff = P_CC * payoffs_A[0] + P_CD * payoffs_A[1] + P_DC * payoffs_A[2] + P_DD * payoffs_A[3]

    # Print the results and the final equation
    print("When both players choose the quantum strategy 'Q', the classical dilemma is resolved.")
    print(f"The probabilities of the outcomes (CC, CD, DC, DD) are ({P_CC:.2f}, {P_CD:.2f}, {P_DC:.2f}, {P_DD:.2f}).")
    print("\nThe payoff for a player at this equilibrium is calculated as:")
    print(f"Payoff = P_CC * {R} + P_CD * {S} + P_DC * {T} + P_DD * {P}")
    print(f"Payoff = {P_CC:.2f} * {R} + {P_CD:.2f} * {S} + {P_DC:.2f} * {T} + {P_DD:.2f} * {P} = {player_A_payoff:.2f}")

    # Return the final numerical answer
    print(f"\n<<<>>>") # Helper for automated result extraction
    print(f"<<<{player_A_payoff}>>>")

run_quantum_prisoners_dilemma()