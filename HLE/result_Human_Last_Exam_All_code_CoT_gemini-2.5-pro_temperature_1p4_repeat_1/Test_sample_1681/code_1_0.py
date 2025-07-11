import numpy as np

def solve_quantum_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol.
    """
    # 1. Define the payoff matrix values
    # (Cooperate, Cooperate)
    R = 5
    # (Cooperate, Defect) for player 1
    S = 0
    # (Defect, Cooperate) for player 1
    T = 7
    # (Defect, Defect)
    P = 1

    # 2. Define quantum states and operators
    # Basis vectors |0> (Cooperate) and |1> (Defect)
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)
    
    # Initial state |00>
    psi_initial = np.kron(ket0, ket0)

    # Standard operators
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    # Entangling operator J for maximal entanglement
    J = (1/np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T

    # Define the quantum strategy 'Q'
    # Q = i * sigma_z
    Q = np.array([[1j, 0], [0, -1j]], dtype=complex)
    
    # Both players choose the quantum strategy Q
    U_A = Q
    U_B = Q

    # 3. Calculate the final state of the game
    # Apply operators: J† * (U_A ⊗ U_B) * J
    final_operator = J_dag @ np.kron(U_A, U_B) @ J
    psi_final = final_operator @ psi_initial

    # 4. Calculate the probability of each classical outcome
    # P_cc = |<00|ψ_final>|^2
    p_cc = np.abs(np.vdot(np.kron(ket0, ket0), psi_final))**2
    # P_cd = |<01|ψ_final>|^2
    p_cd = np.abs(np.vdot(np.kron(ket0, ket1), psi_final))**2
    # P_dc = |<10|ψ_final>|^2
    p_dc = np.abs(np.vdot(np.kron(ket1, ket0), psi_final))**2
    # P_dd = |<11|ψ_final>|^2
    p_dd = np.abs(np.vdot(np.kron(ket1, ket1), psi_final))**2
    
    # 5. Calculate the expected payoff for Player A
    payoff_A = p_cc * R + p_cd * S + p_dc * T + p_dd * P

    # 6. Print the results and the final equation
    print("When both players choose the quantum strategy Q, the probabilities of the outcomes are:")
    print(f"P(Cooperate, Cooperate) = {p_cc:.2f}")
    print(f"P(Cooperate, Defect)   = {p_cd:.2f}")
    print(f"P(Defect,   Cooperate) = {p_dc:.2f}")
    print(f"P(Defect,   Defect)   = {p_dd:.2f}")
    print("\nThis means the outcome is certain cooperation from both sides.")
    
    print("\nThe equilibrium payoff for each player is calculated as follows:")
    # Print the full equation with the numbers
    print(f"Payoff = P(CC)*{R} + P(CD)*{S} + P(DC)*{T} + P(DD)*{P}")
    print(f"Payoff = {p_cc:.2f}*{R} + {p_cd:.2f}*{S} + {p_dc:.2f}*{T} + {p_dd:.2f}*{P}")
    print(f"Payoff = {payoff_A:.2f}")

    print("\nThus, the equilibrium point yields a payoff of 5 for each player, resolving the dilemma.")
    
solve_quantum_dilemma()