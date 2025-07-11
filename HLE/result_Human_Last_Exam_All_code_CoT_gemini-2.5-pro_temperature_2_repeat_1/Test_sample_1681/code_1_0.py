import numpy as np

def calculate_equilibrium():
    """
    Calculates the equilibrium point for the quantum prisoner's dilemma
    using the Eisert-Wilkens-Lewenstein protocol.
    """
    # Payoff Matrix Values
    # R: Reward, S: Sucker, T: Temptation, P: Punishment
    R = 5  # Payoff for (C, C)
    S_A = 0  # Alice's payoff for (C, D)
    T_A = 7  # Alice's payoff for (D, C)
    S_B = 7  # Bob's payoff for (C, D)
    T_B = 0  # Bob's payoff for (D, C)
    P = 1  # Payoff for (D, D)

    # --- Quantum Game Setup ---

    # Define basis states |0> (Cooperate) and |1> (Defect)
    state_0 = np.array([[1], [0]])
    state_1 = np.array([[0], [1]])

    # Initial state of the game |00>
    psi_00 = np.kron(state_0, state_0)

    # Define strategy operators
    # C (Cooperate) is the Identity matrix
    C_op = np.identity(2)
    # D (Defect) is the Pauli-X matrix
    D_op = np.array([[0, 1], [1, 0]])
    # Q (Quantum) is a phase shift i*sigma_z
    Q_op = np.array([[1j, 0], [0, -1j]])

    # Entangling operator J for maximal entanglement
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))

    # Entangled initial state
    psi_initial = J @ psi_00

    # --- Analyze the Nash Equilibrium (D, Q) ---

    # Player strategies: Alice plays Defect (D), Bob plays Quantum (Q)
    U_A = D_op
    U_B = Q_op
    
    # Combined strategy operator
    U_combined = np.kron(U_A, U_B)
    
    # Final state of the game
    psi_final = U_combined @ psi_initial

    # Define basis states for measurement projection
    basis_00 = np.kron(state_0, state_0).T.conj() # <00|
    basis_01 = np.kron(state_0, state_1).T.conj() # <01|
    basis_10 = np.kron(state_1, state_0).T.conj() # <10|
    basis_11 = np.kron(state_1, state_1).T.conj() # <11|
    
    # Calculate probabilities of each classical outcome
    p_cc = np.abs(basis_00 @ psi_final)**2
    p_cd = np.abs(basis_01 @ psi_final)**2
    p_dc = np.abs(basis_10 @ psi_final)**2
    p_dd = np.abs(basis_11 @ psi_final)**2
    
    # Get the single value from the 1x1 matrix result
    p_cc_val, p_cd_val, p_dc_val, p_dd_val = p_cc[0,0], p_cd[0,0], p_dc[0,0], p_dd[0,0]

    # Calculate expected payoffs
    payoff_A = R * p_cc_val + S_A * p_cd_val + T_A * p_dc_val + P * p_dd_val
    payoff_B = R * p_cc_val + S_B * p_cd_val + T_B * p_dc_val + P * p_dd_val

    print("For this payoff matrix, a Nash Equilibrium exists at the strategy pair (Defect, Quantum).")
    print("The final state probabilities are:")
    print(f"P(C,C) = {p_cc_val:.2f}")
    print(f"P(C,D) = {p_cd_val:.2f}")
    print(f"P(D,C) = {p_dc_val:.2f}")
    print(f"P(D,D) = {p_dd_val:.2f}\n")
    
    print("Alice's Payoff (Player 1 playing D):")
    print(f"Payoff_A = R*P(C,C) + S_A*P(C,D) + T_A*P(D,C) + P*P(D,D)")
    print(f"Payoff_A = {R} * {p_cc_val:.2f} + {S_A} * {p_cd_val:.2f} + {T_A} * {p_dc_val:.2f} + {P} * {p_dd_val:.2f} = {payoff_A:.1f}\n")

    print("Bob's Payoff (Player 2 playing Q):")
    print(f"Payoff_B = R*P(C,C) + S_B*P(C,D) + T_B*P(D,C) + P*P(D,D)")
    print(f"Payoff_B = {R} * {p_cc_val:.2f} + {S_B} * {p_cd_val:.2f} + {T_B} * {p_dc_val:.2f} + {P} * {p_dd_val:.2f} = {payoff_B:.1f}\n")

    print(f"The equilibrium point results in a payoff of {payoff_A:.1f} for each player.")

calculate_equilibrium()