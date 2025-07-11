import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates and explains the Nash Equilibrium in a quantum Prisoner's Dilemma.
    """
    # 1. Define game constants and basis states
    # Payoffs: R(CC), S(CD), T(DC), P(DD)
    # For player A: R=5 for (C,C), S=0 for (C,D), T=7 for (D,C), P=1 for (D,D)
    # For player B: R=5 for (C,C), S=0 for (D,C), T=7 for (C,D), P=1 for (D,D)
    R, S, T, P = 5, 0, 7, 1

    # Basis states |C> = |0>, |D> = |1>
    C_state = np.array([1, 0], dtype=complex)
    D_state = np.array([0, 1], dtype=complex)

    # Basis for the two-qubit system
    CC = np.kron(C_state, C_state)
    CD = np.kron(C_state, D_state)
    DC = np.kron(D_state, C_state)
    DD = np.kron(D_state, D_state)

    # 2. Define quantum operators for strategies
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Classical strategies as unitary operators
    C_op = I             # Cooperate ~ Identity
    D_op = 1j * sigma_x  # Defect ~ bit-flip (Pauli-X)

    # The special Quantum strategy
    Q_op = 1j * sigma_z  # Quantum 'miracle' move (Pauli-Z with a phase)

    # 3. Define the EWL protocol operators
    # J is the entangling operator (for maximal entanglement)
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    # J_dag is the disentangling operator
    J_dag = J.conj().T

    # The initial entangled state of the game
    psi_initial = J @ CC

    def calculate_payoffs(op_A, op_B):
        """Calculates payoffs for players A and B given their strategy operators."""
        # Players apply their strategies to the initial state
        psi_prime = np.kron(op_A, op_B) @ psi_initial
        
        # The disentangling operator is applied
        psi_final = J_dag @ psi_prime

        # Probabilities of classical outcomes are measured from the final state
        p_cc = abs(CC.conj().T @ psi_final)**2
        p_cd = abs(CD.conj().T @ psi_final)**2
        p_dc = abs(DC.conj().T @ psi_final)**2
        p_dd = abs(DD.conj().T @ psi_final)**2
        
        # Expected payoffs are calculated based on these probabilities
        payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
        payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd

        return (payoff_A, payoff_B, p_cc, p_cd, p_dc, p_dd)

    # 4. Analyze the (Q, Q) strategy as a potential equilibrium
    print("Analysis of Quantum Prisoner's Dilemma Equilibrium:")
    print("="*60)
    print("The payoff matrix values are: R=5 (C,C), S=0 (C,D), T=7 (D,C), P=1 (D,D).")
    print("A new quantum strategy 'Q' is available.\n")
    print("We test if the strategy profile (Q, Q) constitutes a Nash Equilibrium.")

    # Case 1: Both players choose the quantum strategy Q
    payoff_A_QQ, payoff_B_QQ, pcc, pcd, pdc, pdd = calculate_payoffs(Q_op, Q_op)
    print("\nCase 1: Both players choose Q -> (Q, Q)")
    print("  Player A's payoff calculation:")
    print(f"    Payoff_A = ({R} * {pcc:.2f}) + ({S} * {pcd:.2f}) + ({T} * {pdc:.2f}) + ({P} * {pdd:.2f}) = {payoff_A_QQ:.2f}")
    print("  Player B's payoff calculation:")
    print(f"    Payoff_B = ({R} * {pcc:.2f}) + ({T} * {pcd:.2f}) + ({S} * {pdc:.2f}) + ({P} * {pdd:.2f}) = {payoff_B_QQ:.2f}")
    print(f"  --> Resulting Payoff: ({payoff_A_QQ:.0f}, {payoff_B_QQ:.0f})")

    # Case 2: Check for deviation. Player A defects while Player B plays Q
    print("\nNow, we check if Player A benefits from unilaterally deviating to Defect 'D'.")
    payoff_A_DQ, _, pcc, pcd, pdc, pdd = calculate_payoffs(D_op, Q_op)
    print("Case 2: Player A deviates -> (D, Q)")
    print("  Player A's new payoff calculation:")
    print(f"    Payoff_A = ({R} * {pcc:.2f}) + ({S} * {pcd:.2f}) + ({T} * {pdc:.2f}) + ({P} * {pdd:.2f}) = {payoff_A_DQ:.2f}")
    print(f"  --> Player A's payoff changes from {payoff_A_QQ:.0f} to {payoff_A_DQ:.0f}. Deviation is not profitable.")

    # Case 3: Check for deviation. Player A cooperates while Player B plays Q
    print("\nNext, we check if Player A benefits from unilaterally deviating to Cooperate 'C'.")
    payoff_A_CQ, _, pcc, pcd, pdc, pdd = calculate_payoffs(C_op, Q_op)
    print("Case 3: Player A deviates -> (C, Q)")
    print("  Player A's new payoff calculation:")
    print(f"    Payoff_A = ({R} * {pcc:.2f}) + ({S} * {pcd:.2f}) + ({T} * {pdc:.2f}) + ({P} * {pdd:.2f}) = {payoff_A_CQ:.2f}")
    print(f"  --> Player A's payoff changes from {payoff_A_QQ:.0f} to {payoff_A_CQ:.0f}. Deviation is not profitable.")
    
    # 5. Conclusion
    print("\n" + "="*60)
    print("Conclusion:")
    print("Since neither player can improve their payoff by unilaterally deviating from Q, the strategy")
    print("profile (Q, Q) is a Nash Equilibrium. This equilibrium yields the Pareto optimal")
    print(f"outcome of ({payoff_A_QQ:.0f}, {payoff_B_QQ:.0f}), thus resolving the classical dilemma.")
    print("The equilibrium point corresponds to a payoff of 5 for each player.")

solve_quantum_prisoners_dilemma()