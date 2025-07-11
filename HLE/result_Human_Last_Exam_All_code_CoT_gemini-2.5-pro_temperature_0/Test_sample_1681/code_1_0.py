import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    with maximal entanglement.
    """
    # Payoff matrix values:
    # R: Reward (C,C), S: Sucker (C,D), T: Temptation (D,C), P: Punishment (D,D)
    R, S, T, P = 5, 0, 7, 1

    # --- Quantum Game Setup ---

    # Identity and Pauli-X gates
    I = np.eye(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # Entangling operator J for maximal entanglement (gamma = pi/2)
    gamma = np.pi / 2
    J = (1/np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dagger = J.conj().T

    # Initial state is |00>
    psi_00 = np.array([1, 0, 0, 0], dtype=complex)

    # --- Equilibrium Strategies ---

    # The new equilibrium strategy 'Q' (a specific unitary operator)
    # Q = i * sigma_z
    Q = np.array([[1j, 0], [0, -1j]], dtype=complex)

    # Both players adopt the optimal quantum strategy Q
    U_A = Q
    U_B = Q

    # The combined strategy operator for both players
    U_total = np.kron(U_A, U_B)

    # --- Calculation ---

    # 1. The initial state |00> is entangled by J
    # 2. Players apply their strategies U_A and U_B
    # 3. The state is disentangled by J_dagger before measurement
    psi_final = J_dagger @ U_total @ J @ psi_00

    # The probabilities of the classical outcomes are the squared magnitudes
    # of the amplitudes of the final state vector.
    # The basis is |00>, |01>, |10>, |11>
    probabilities = np.abs(psi_final)**2
    P_CC, P_CD, P_DC, P_DD = probabilities[0], probabilities[1], probabilities[2], probabilities[3]

    # Calculate the expected payoff for each player
    payoff_A = P_CC * R + P_CD * S + P_DC * T + P_DD * P
    payoff_B = P_CC * R + P_CD * T + P_DC * S + P_DD * P

    # --- Output the Result ---
    print("At the quantum equilibrium point:")
    print(f"The probability of (Cooperate, Cooperate) is {P_CC:.2f}")
    print(f"The probability of (Cooperate, Defect) is {P_CD:.2f}")
    print(f"The probability of (Defect, Cooperate) is {P_DC:.2f}")
    print(f"The probability of (Defect, Defect) is {P_DD:.2f}")
    print("\nThe final payoff calculation is:")
    
    # Print the full equation for each player's payoff
    print(f"Player A Payoff = ({P_CC:.2f} * {R}) + ({P_CD:.2f} * {S}) + ({P_DC:.2f} * {T}) + ({P_DD:.2f} * {P}) = {payoff_A:.1f}")
    print(f"Player B Payoff = ({P_CC:.2f} * {R}) + ({P_CD:.2f} * {T}) + ({P_DC:.2f} * {S}) + ({P_DD:.2f} * {P}) = {payoff_B:.1f}")
    
    print(f"\nThe equilibrium point is a payoff of ({payoff_A:.1f}, {payoff_B:.1f}).")


solve_quantum_prisoners_dilemma()