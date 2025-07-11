import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum prisoner's dilemma.
    """
    # Payoff matrix values: (R, S, T, P)
    # (C,C) -> (R,R) = (5,5)
    # (C,D) -> (S,T) = (0,7)
    # (D,C) -> (T,S) = (7,0)
    # (D,D) -> (P,P) = (1,1)
    R, S, T, P = 5, 0, 7, 1

    # Define Pauli-X and Identity matrices
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # Entangling operator J for maximal entanglement (gamma = pi/2)
    J = (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x)) / np.sqrt(2)
    J_H = J.conj().T # Hermitian conjugate

    # Define the strategy operator U(theta, phi)
    def U(theta, phi):
        return np.array([
            [np.cos(theta / 2), 1j * np.exp(-1j * phi) * np.sin(theta / 2)],
            [1j * np.exp(1j * phi) * np.sin(theta / 2), np.cos(theta / 2)]
        ], dtype=complex)

    # The Nash Equilibrium strategy in the quantum game is Q = U(pi/2, 0)
    theta_eq = np.pi / 2
    phi_eq = 0.0

    U_A = U(theta_eq, phi_eq)
    U_B = U(theta_eq, phi_eq)

    # Initial state is |CC> or |00>
    psi_initial = np.array([1, 0, 0, 0], dtype=complex)

    # The full evolution operator
    # Note: The standard protocol applies J, then strategies, then J_H
    full_evolution = J_H @ np.kron(U_A, U_B) @ J

    # Calculate the final state
    psi_final = full_evolution @ psi_initial

    # Basis states for measurement
    C00 = np.array([1, 0, 0, 0], dtype=complex) # Cooperate, Cooperate
    C01 = np.array([0, 1, 0, 0], dtype=complex) # Cooperate, Defect
    C10 = np.array([0, 0, 1, 0], dtype=complex) # Defect, Cooperate
    C11 = np.array([0, 0, 0, 1], dtype=complex) # Defect, Defect

    # Calculate probabilities of each outcome
    P_cc = np.abs(np.vdot(C00, psi_final))**2
    P_cd = np.abs(np.vdot(C01, psi_final))**2
    P_dc = np.abs(np.vdot(C10, psi_final))**2
    P_dd = np.abs(np.vdot(C11, psi_final))**2

    # Calculate the payoffs for each player
    payoff_A = R * P_cc + S * P_cd + T * P_dc + P * P_dd
    payoff_B = R * P_cc + T * P_cd + S * P_dc + P * P_dd

    print("Quantum Prisoner's Dilemma Equilibrium Calculation")
    print("================================================")
    print(f"The Nash Equilibrium strategy is Q = U(theta=pi/2, phi=0).")
    print(f"When both players choose Q, the outcome probabilities are:")
    print(f"P(Cooperate, Cooperate) = {P_cc:.4f}")
    print(f"P(Cooperate, Defect)   = {P_cd:.4f}")
    print(f"P(Defect, Cooperate)   = {P_dc:.4f}")
    print(f"P(Defect, Defect)     = {P_dd:.4f}")
    print("\nPlayer A's equilibrium payoff is calculated as:")
    print(f"Payoff_A = (R * P_cc) + (S * P_cd) + (T * P_dc) + (P * P_dd)")
    print(f"Payoff_A = ({R} * {P_cc:.2f}) + ({S} * {P_cd:.2f}) + ({T} * {P_dc:.2f}) + ({P} * {P_dd:.2f}) = {payoff_A:.2f}")
    
    print("\nPlayer B's equilibrium payoff is calculated as:")
    print(f"Payoff_B = (R * P_cc) + (T * P_cd) + (S * P_dc) + (P * P_dd)")
    print(f"Payoff_B = ({R} * {P_cc:.2f}) + ({T} * {P_cd:.2f}) + ({S} * {P_dc:.2f}) + ({P} * {P_dd:.2f}) = {payoff_B:.2f}")
    
    print(f"\nThe equilibrium point results in a payoff of {payoff_A:.2f} for each player.")

solve_quantum_prisoners_dilemma()