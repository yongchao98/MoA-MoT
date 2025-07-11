import numpy as np

def run_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma.
    """
    # Payoff matrix [ [(A_CC, B_CC), (A_CD, B_CD)], [(A_DC, B_DC), (A_DD, B_DD)] ]
    payoff_A = np.array([[5, 0], [7, 1]])
    payoff_B = np.array([[5, 7], [0, 1]])

    # Define standard matrices and quantum states
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)

    ket_0 = np.array([[1], [0]], dtype=complex)
    ket_1 = np.array([[0], [1]], dtype=complex)
    ket_00 = np.kron(ket_0, ket_0) # Represents (Cooperate, Cooperate) initial state

    # Step 1: Define the entangling operator J
    # J = exp(i * gamma/2 * (sigma_x ⊗ sigma_x))
    # We choose an optimal scenario with maximum entanglement, gamma = pi/2.
    gamma = np.pi / 2
    J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(sigma_x, sigma_x)
    J_dagger = J.conj().T

    # Step 2: Define player strategies U(theta)
    # U(theta) = exp(i * theta/2 * sigma_y)
    def U(theta):
        return np.cos(theta / 2) * I + 1j * np.sin(theta / 2) * sigma_y

    # The new Nash Equilibrium in the quantum game is the strategy theta = pi/2 ("Q")
    # C (Cooperate) is U(0) and D (Defect) is U(pi).
    theta_A = np.pi / 2 # Player A chooses strategy Q
    theta_B = np.pi / 2 # Player B chooses strategy Q

    U_A = U(theta_A)
    U_B = U(theta_B)

    # Step 3: Calculate the final state
    # psi_final = J_dagger * (U_A ⊗ U_B) * J * |00>
    U_kron = np.kron(U_A, U_B)
    psi_final = J_dagger @ U_kron @ J @ ket_00

    # Step 4: Calculate the probabilities of each outcome
    # P_ij = |<ij|psi_final>|^2
    prob_00 = np.abs(np.kron(ket_0, ket_0).conj().T @ psi_final)**2
    prob_01 = np.abs(np.kron(ket_0, ket_1).conj().T @ psi_final)**2
    prob_10 = np.abs(np.kron(ket_1, ket_0).conj().T @ psi_final)**2
    prob_11 = np.abs(np.kron(ket_1, ket_1).conj().T @ psi_final)**2

    P_CC = prob_00[0,0]
    P_CD = prob_01[0,0]
    P_DC = prob_10[0,0]
    P_DD = prob_11[0,0]

    # Step 5: Calculate the expected payoffs
    PA_CC, PA_CD, PA_DC, PA_DD = payoff_A.flatten()
    PB_CC, PB_CD, PB_DC, PB_DD = payoff_B.flatten()
    
    payoff_A_eq = PA_CC * P_CC + PA_CD * P_CD + PA_DC * P_DC + PA_DD * P_DD
    payoff_B_eq = PB_CC * P_CC + PB_CD * P_CD + PB_DC * P_DC + PB_DD * P_DD

    print("At the (Q, Q) equilibrium:")
    print("-" * 30)
    print(f"Probabilities (CC, CD, DC, DD): ({P_CC:.3f}, {P_CD:.3f}, {P_DC:.3f}, {P_DD:.3f})")
    print("\nPayoff Calculation for Player A:")
    print(f"E_A = ({PA_CC} * {P_CC:.3f}) + ({PA_CD} * {P_CD:.3f}) + ({PA_DC} * {P_DC:.3f}) + ({PA_DD} * {P_DD:.3f})")
    print(f"E_A = {payoff_A_eq:.3f}")
    
    print("\nPayoff Calculation for Player B:")
    print(f"E_B = ({PB_CC} * {P_CC:.3f}) + ({PB_CD} * {P_CD:.3f}) + ({PB_DC} * {P_DC:.3f}) + ({PB_DD} * {P_DD:.3f})")
    print(f"E_B = {payoff_B_eq:.3f}")
    
    print("\nThe equilibrium point payoff is ({:.1f}, {:.1f}).".format(payoff_A_eq, payoff_B_eq))


run_quantum_prisoners_dilemma()