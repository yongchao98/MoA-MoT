import numpy as np
from scipy.linalg import expm

def solve_quantum_pd():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol.
    """

    # 1. Define the payoff matrix values
    # Cooperate-Cooperate (Reward)
    R = 5
    # Cooperate-Defect (Sucker)
    S = 0
    # Defect-Cooperate (Temptation)
    T = 7
    # Defect-Defect (Punishment)
    P = 1

    # 2. Define the quantum operators
    # Pauli-X matrix
    sigma_x = np.array([[0, 1], [1, 0]])

    # Entangling operator J for maximal entanglement (gamma = pi/2)
    # J = exp(-i * (pi/4) * (sigma_x ⊗ sigma_x))
    gamma = np.pi / 2
    J = expm(-1j * (gamma / 2) * np.kron(sigma_x, sigma_x))
    
    # The disentangling operator is the inverse (conjugate transpose) of J
    J_dag = J.conj().T

    # 3. Define the players' strategies (unitary operators)
    # The new quantum Nash Equilibrium strategy, Q
    # This corresponds to U(theta=pi/2, phi=pi/2) in the Eisert et al. parametrization.
    # U(theta, phi) = [[e^(i*phi)*cos(theta/2), sin(theta/2)],
    #                  [-sin(theta/2), e^(-i*phi)*cos(theta/2)]]
    theta_q = np.pi / 2
    phi_q = np.pi / 2
    c = np.cos(theta_q / 2)
    s = np.sin(theta_q / 2)
    Q = np.array([
        [np.exp(1j * phi_q) * c, s],
        [-s, np.exp(-1j * phi_q) * c]
    ])

    # In the equilibrium, both players choose strategy Q.
    U_A = Q
    U_B = Q

    # 4. Calculate the evolution of the quantum state
    # Initial state is |00>
    psi_initial = np.array([1, 0, 0, 0])

    # The full evolution operator for the (Q, Q) strategy pair
    evolution_op = J_dag @ np.kron(U_A, U_B) @ J
    
    # The final state of the system
    psi_final = evolution_op @ psi_initial

    # 5. Calculate the probabilities of the four classical outcomes
    # P_ij = |<ij|psi_final>|^2
    # |00> -> CC, |01> -> CD, |10> -> DC, |11> -> DD
    P_CC = np.abs(psi_final[0])**2
    P_CD = np.abs(psi_final[1])**2
    P_DC = np.abs(psi_final[2])**2
    P_DD = np.abs(psi_final[3])**2

    # 6. Calculate the expected payoffs for each player
    payoff_A = P_CC * R + P_CD * S + P_DC * T + P_DD * P
    payoff_B = P_CC * R + P_CD * T + P_DC * S + P_DD * P

    print("Quantum Prisoner's Dilemma Equilibrium Calculation:")
    print("-" * 50)
    print(f"Payoff Matrix: R={R}, S={S}, T={T}, P={P}")
    print(f"Equilibrium Strategy: Both players choose the quantum strategy 'Q'.")
    print("\nResulting Probabilities:")
    print(f"P(Cooperate, Cooperate) = {P_CC:.4f}")
    print(f"P(Cooperate, Defect)   = {P_CD:.4f}")
    print(f"P(Defect, Cooperate)   = {P_DC:.4f}")
    print(f"P(Defect, Defect)      = {P_DD:.4f}")

    # For the (Q, Q) equilibrium, payoffs are symmetric.
    print("\nEquilibrium Payoff Equation for Player A:")
    # Using 'round' to handle potential floating point inaccuracies for clean output
    print(f"{payoff_A:.2f} = {round(P_CC,2)} * {R} + {round(P_CD,2)} * {S} + {round(P_DC,2)} * {T} + {round(P_DD,2)} * {P}")


solve_quantum_pd()