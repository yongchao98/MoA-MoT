import numpy as np
from scipy.linalg import expm

def solve_quantum_prisoners_dilemma():
    """
    Solves the Prisoner's Dilemma in a quantum setting using the EWL protocol.
    
    This function calculates the payoffs for the Nash Equilibrium that resolves
    the classical dilemma.
    """

    # Payoff matrix values
    # R: Reward (C,C), S: Sucker (C,D), T: Temptation (D,C), P: Punishment (D,D)
    R = 5
    S = 0
    T = 7
    P = 1
    
    payoffs = {'R': R, 'S': S, 'T': T, 'P': P}

    # Define Pauli matrices and Identity
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # 1. Initial State Setup
    # Create the entangling operator J for maximal entanglement (gamma = pi/2)
    # J = expm(i * (pi/4) * (sigma_y kronecker sigma_y))
    J_exponent = 1j * (np.pi / 4) * np.kron(sigma_y, sigma_y)
    J = expm(J_exponent)
    
    # The initial classical state is |00> (Cooperate, Cooperate)
    psi_initial = np.array([1, 0, 0, 0], dtype=complex)

    # 2. Player Strategies
    # C = Cooperate (Identity)
    # Q = The new quantum strategy that forms the equilibrium
    C = I
    Q = 1j * sigma_z
    
    # We analyze the equilibrium case where both players choose Q
    U_A = Q
    U_B = Q

    # 3. Final State Calculation
    # Apply the full protocol: J_dagger * (U_A x U_B) * J * psi_initial
    J_dagger = J.conj().T
    U_total = np.kron(U_A, U_B)
    
    psi_final = J_dagger @ U_total @ J @ psi_initial
    
    # 4. Calculate Probabilities and Payoffs
    # The probabilities of measuring the four basis states |00>, |01>, |10>, |11>
    probs = np.abs(psi_final)**2
    
    p_cc = probs[0] # P(|00>) - Cooperate, Cooperate
    p_cd = probs[1] # P(|01>) - Cooperate, Defect
    p_dc = probs[2] # P(|10>) - Defect, Cooperate
    p_dd = probs[3] # P(|11>) - Defect, Defect

    # Calculate expected payoffs for Alice and Bob
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd

    # Output the results
    print("Quantum Prisoner's Dilemma Equilibrium Calculation")
    print("==================================================")
    print("Chosen Equilibrium Strategy Pair: (Q, Q) where Q = i * sigma_z\n")
    print(f"Probabilities of outcomes:")
    print(f"P(CC) = {p_cc:.2f}, P(CD) = {p_cd:.2f}, P(DC) = {p_dc:.2f}, P(DD) = {p_dd:.2f}\n")
    print("Player A's payoff calculation:")
    print(f"{R} * {p_cc:.2f} + {S} * {p_cd:.2f} + {T} * {p_dc:.2f} + {P} * {p_dd:.2f} = {payoff_A:.2f}")

    print("\nPlayer B's payoff calculation:")
    print(f"{R} * {p_cc:.2f} + {T} * {p_cd:.2f} + {S} * {p_dc:.2f} + {P} * {p_dd:.2f} = {payoff_B:.2f}")

    print("\nThus, the equilibrium point results in a payoff of (5.00, 5.00), resolving the dilemma.")


solve_quantum_prisoners_dilemma()