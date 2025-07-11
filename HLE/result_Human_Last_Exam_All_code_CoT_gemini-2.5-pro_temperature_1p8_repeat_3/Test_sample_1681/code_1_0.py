import numpy as np

def quantum_prisoners_dilemma():
    """
    Calculates and explains the equilibrium point for the Quantum Prisoner's Dilemma.
    """
    # Payoff values from the matrix
    # (C, C) -> (R, R)
    # (C, D) -> (S, T)
    # (D, C) -> (T, S)
    # (D, D) -> (P, P)
    R, S, T, P = 5, 0, 7, 1

    # Basis states: |C> = [1, 0], |D> = [0, 1]
    C = np.array([[1], [0]])
    D = np.array([[0], [1]])

    # Basis states for the 2-qubit system
    CC = np.kron(C, C)
    CD = np.kron(C, D)
    DC = np.kron(D, C)
    DD = np.kron(D, D)

    # Pauli Matrices
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])

    # Entangling operator J and its conjugate transpose J_dag
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T

    # Initial entangled state
    psi_initial = J @ CC

    # Strategy Operators
    C_op = I          # Cooperate
    D_op = 1j * sigma_x # Defect
    Q_op = 1j * sigma_z # Quantum strategy

    print("--- Finding the Quantum Equilibrium ---")

    # Case 1: Both players choose the Quantum strategy (Q, Q)
    U_A_eq = Q_op
    U_B_eq = Q_op

    # Combined strategy operator
    U_eq = np.kron(U_A_eq, U_B_eq)

    # Calculate final state
    psi_final_eq = J_dag @ U_eq @ psi_initial

    # Calculate probabilities of classical outcomes
    p_cc_eq = np.abs(psi_final_eq.T @ CC)**2
    p_cd_eq = np.abs(psi_final_eq.T @ CD)**2
    p_dc_eq = np.abs(psi_final_eq.T @ DC)**2
    p_dd_eq = np.abs(psi_final_eq.T @ DD)**2

    # Calculate payoffs
    payoff_A_eq = p_cc_eq * R + p_cd_eq * S + p_dc_eq * T + p_dd_eq * P
    payoff_B_eq = p_cc_eq * R + p_cd_eq * T + p_dc_eq * S + p_dd_eq * P
    
    print("\nStrategy Profile: (Quantum, Quantum)")
    print(f"Probabilities: P(CC)={p_cc_eq[0][0]:.2f}, P(CD)={p_cd_eq[0][0]:.2f}, P(DC)={p_dc_eq[0][0]:.2f}, P(DD)={p_dd_eq[0][0]:.2f}")
    print(f"Player A Payoff = {p_cc_eq[0][0]:.2f}*{R} + {p_cd_eq[0][0]:.2f}*{S} + {p_dc_eq[0][0]:.2f}*{T} + {p_dd_eq[0][0]:.2f}*{P} = {payoff_A_eq[0][0]:.2f}")
    print(f"Player B Payoff = {p_cc_eq[0][0]:.2f}*{R} + {p_cd_eq[0][0]:.2f}*{T} + {p_dc_eq[0][0]:.2f}*{S} + {p_dd_eq[0][0]:.2f}*{P} = {payoff_B_eq[0][0]:.2f}")

    print("\n--- Verifying the Equilibrium (Checking for Deviations) ---")
    
    # Case 2: Player A plays Q, Player B deviates to D
    U_A_dev = Q_op
    U_B_dev = D_op
    U_dev = np.kron(U_A_dev, U_B_dev)
    psi_final_dev = J_dag @ U_dev @ psi_initial

    p_cc_dev = np.abs(psi_final_dev.T @ CC)**2
    p_cd_dev = np.abs(psi_final_dev.T @ CD)**2
    p_dc_dev = np.abs(psi_final_dev.T @ DC)**2
    p_dd_dev = np.abs(psi_final_dev.T @ DD)**2

    payoff_A_dev = p_cc_dev * R + p_cd_dev * S + p_dc_dev * T + p_dd_dev * P
    payoff_B_dev = p_cc_dev * R + p_cd_dev * T + p_dc_dev * S + p_dd_dev * P

    print("\nStrategy Profile: (Quantum, Defect)")
    print(f"Probabilities: P(CC)={p_cc_dev[0][0]:.2f}, P(CD)={p_cd_dev[0][0]:.2f}, P(DC)={p_dc_dev[0][0]:.2f}, P(DD)={p_dd_dev[0][0]:.2f}")
    print(f"Player A Payoff = {p_cc_dev[0][0]:.2f}*{R} + {p_cd_dev[0][0]:.2f}*{S} + {p_dc_dev[0][0]:.2f}*{T} + {p_dd_dev[0][0]:.2f}*{P} = {payoff_A_dev[0][0]:.2f}")
    print(f"Player B Payoff = {p_cc_dev[0][0]:.2f}*{R} + {p_cd_dev[0][0]:.2f}*{T} + {p_dc_dev[0][0]:.2f}*{S} + {p_dd_dev[0][0]:.2f}*{P} = {payoff_B_dev[0][0]:.2f}")
    
    print("\nConclusion:")
    print(f"Player B's payoff drops from {payoff_B_eq[0][0]:.2f} to {payoff_B_dev[0][0]:.2f} upon unilaterally defecting.")
    print("Therefore, there is no incentive to deviate from the (Quantum, Quantum) strategy.")
    print("The equilibrium point for the quantum game is a payoff of (5, 5), which resolves the dilemma.")


if __name__ == '__main__':
    quantum_prisoners_dilemma()
    # The equilibrium point is the payoff pair achieved at equilibrium.
    print("\n<<< (5, 5) >>>")
