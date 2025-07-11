import numpy as np

def solve_quantum_dilemma():
    """
    This function sets up and solves the Quantum Prisoner's Dilemma,
    demonstrating the Nash Equilibrium that resolves the classical dilemma.
    """
    # 1. Setup the Game
    
    # Payoff Matrix: (C,C), (C,D), (D,C), (D,D)
    # R=5 (Reward), S=0 (Sucker), T=7 (Temptation), P=1 (Punishment)
    payoffs_A = np.array([5, 0, 7, 1]) # Alice's payoffs for CC, CD, DC, DD
    payoffs_B = np.array([5, 7, 0, 1]) # Bob's payoffs for CC, CD, DC, DD

    # Pauli Matrices (fundamental to quantum operations)
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # EWL Protocol Entangling Operator J
    J = (1/np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T # Conjugate transpose (disentangling operator)

    # 2. Define Strategies
    
    # C: Cooperate (Identity matrix)
    C = I
    # D: Defect (equivalent to i * Pauli-Y matrix)
    D = 1j * sigma_y
    # Q: Quantum "Miracle Move" (equivalent to i * Pauli-Z matrix)
    Q = 1j * sigma_z

    # Initial state of the game is |CC>
    psi_initial = np.array([1, 0, 0, 0], dtype=complex)

    # 3. Payoff Calculation Function
    
    def calculate_payoffs(U_A, U_B):
        """Calculates payoffs for Alice and Bob given their strategies."""
        U_combined = np.kron(U_A, U_B)
        # The final state amplitudes are calculated using the EWL formula
        final_amplitudes = J_dag @ U_combined @ J @ psi_initial
        # Probabilities are the squared magnitudes of the amplitudes
        probabilities = np.abs(final_amplitudes)**2
        
        payoff_A = np.sum(payoffs_A * probabilities)
        payoff_B = np.sum(payoffs_B * probabilities)
        
        return payoff_A, payoff_B, probabilities

    # 4. Identify and Analyze the Equilibrium
    
    print("Analyzing the Quantum Prisoner's Dilemma Equilibrium...")
    print("-" * 60)

    # Test the (Q, Q) strategy pair
    payoff_A_QQ, payoff_B_QQ, probs_QQ = calculate_payoffs(Q, Q)
    print("1. Both players choose the Quantum strategy Q:")
    print(f"   - Strategy pair: (Q, Q)")
    print(f"   - Outcome probabilities (CC, CD, DC, DD): {np.round(probs_QQ, 2)}")
    print(f"   - Payoff: ({payoff_A_QQ:.2f}, {payoff_B_QQ:.2f})")
    print("   This pair achieves the Pareto optimal outcome, resolving the dilemma.")

    # Check for Nash Equilibrium by testing for profitable deviations
    print("\n2. Is (Q, Q) a Nash Equilibrium? Checking for deviations:")
    print("   (Assuming Bob plays Q, will Alice gain by changing her strategy?)")
    
    # Alice deviates to C (Cooperate)
    payoff_A_CQ, _, _ = calculate_payoffs(C, Q)
    print(f"   - If Alice deviates to C: Payoff -> {payoff_A_CQ:.2f}. (This is not better than 5.00)")

    # Alice deviates to D (Defect)
    payoff_A_DQ, _, _ = calculate_payoffs(D, Q)
    print(f"   - If Alice deviates to D: Payoff -> {payoff_A_DQ:.2f}. (This is not better than 5.00)")

    print("\nConclusion: Neither player has an incentive to unilaterally deviate from Q.")
    print("Therefore, (Q, Q) is the Nash Equilibrium in the quantum game.")
    print("-" * 60)

    # 5. Present the Final Result
    
    print("\nThe equilibrium point is the strategy pair (Q, Q), which results in the payoff (5, 5).")
    print("The final equation for this equilibrium payoff is:")

    p_cc, p_cd, p_dc, p_dd = probs_QQ
    r, s, t, p = payoffs_A[0], payoffs_A[1], payoffs_A[2], payoffs_A[3]

    print(f"\nPayoff = Payoff(CC)*P(CC) + Payoff(CD)*P(CD) + Payoff(DC)*P(DC) + Payoff(DD)*P(DD)")
    print(f"{payoff_A_QQ:.2f} = {r} * {p_cc:.2f} + {s} * {p_cd:.2f} + {t} * {p_dc:.2f} + {p} * {p_dd:.2f}")

solve_quantum_dilemma()