import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Solves the Prisoner's Dilemma in a quantum setting to find the equilibrium point.
    """
    # Define basis states and Pauli matrices (quantum gates)
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    ket00 = np.array([[1], [0], [0], [0]], dtype=complex)

    # Define the players' strategies as unitary operators
    # Cooperate strategy
    C = I
    # Defect strategy
    D = 1j * sx
    # Quantum "miracle move" strategy
    Q = 1j * sz

    # Setup the game using the EWL protocol
    # J is the entangling operator that creates the initial state from |00>
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sx, sx))
    # J_dag is the disentangling operator
    J_dag = J.conj().T

    # Payoff matrix V = (Alice's payoff, Bob's payoff)
    # The four outcomes are (Cooperate, Cooperate), (Cooperate, Defect),
    # (Defect, Cooperate), and (Defect, Defect).
    payoff_values = {
        'A': np.array([5, 0, 7, 1]), # Payoffs for Alice
        'B': np.array([5, 7, 0, 1])  # Payoffs for Bob
    }

    def calculate_payoffs(U_A, U_B, U_A_name, U_B_name):
        """
        Calculates and returns the expected payoffs for Alice and Bob for given strategies.
        
        The final state is calculated as: |psi_f> = J_dag * (U_A kron U_B) * J * |00>
        """
        U = np.kron(U_A, U_B)
        psi_final = J_dag @ U @ J @ ket00

        # Probabilities are the squared magnitudes of the final state's amplitudes
        probs = np.abs(psi_final.flatten())**2

        # Expected payoffs are the sum of (probability * payoff) for each outcome
        payoff_A = np.sum(probs * payoff_values['A'])
        payoff_B = np.sum(probs * payoff_values['B'])

        return payoff_A, payoff_B, probs

    # --- Analysis ---
    print("Analyzing payoffs for different strategy pairs...")

    # 1. Equilibrium Point: Both players use the quantum strategy (Q, Q)
    payoff_A_QQ, payoff_B_QQ, probs_QQ = calculate_payoffs(Q, Q, 'Q', 'Q')
    print(f"\n1. Payoff(Alice=Q, Bob=Q): ({payoff_A_QQ:.2f}, {payoff_B_QQ:.2f})")
    print("   When both players choose the quantum strategy Q, the outcome is a certain 'Cooperate, Cooperate'.")
    print("   This yields the Pareto-optimal payoff, resolving the dilemma.")


    # 2. Deviation Check 1: Alice deviates to Cooperate (C)
    payoff_A_CQ, payoff_B_CQ, _ = calculate_payoffs(C, Q, 'C', 'Q')
    print(f"\n2. Payoff(Alice=C, Bob=Q): ({payoff_A_CQ:.2f}, {payoff_B_CQ:.2f})")
    print(f"   If Alice deviates to C, her payoff drops from {payoff_A_QQ:.2f} to {payoff_A_CQ:.2f}. No incentive to deviate.")

    # 3. Deviation Check 2: Alice deviates to Defect (D)
    payoff_A_DQ, payoff_B_DQ, _ = calculate_payoffs(D, Q, 'D', 'Q')
    print(f"\n3. Payoff(Alice=D, Bob=Q): ({payoff_A_DQ:.2f}, {payoff_B_DQ:.2f})")
    print(f"   If Alice deviates to D, her payoff drops from {payoff_A_QQ:.2f} to {payoff_A_DQ:.2f}. No incentive to deviate.")
    
    # 4. Classical Nash Equilibrium: (D, D) for reference
    payoff_A_DD, payoff_B_DD, _ = calculate_payoffs(D, D, 'D', 'D')
    print(f"\n4. Classical Payoff(Alice=D, Bob=D): ({payoff_A_DD:.2f}, {payoff_B_DD:.2f})")
    print(f"   The quantum equilibrium payoff of ({payoff_A_QQ:.2f}, {payoff_B_QQ:.2f}) is far superior to the classical Nash Equilibrium.")
    

    print("\n----------------------------------------------------")
    print("EQUILIBRIUM POINT CALCULATION")
    print("----------------------------------------------------")
    print("The equilibrium point is the strategy pair (Q, Q) with the following payoff calculation:")
    
    print("\nAlice's Equilibrium Payoff Calculation:")
    print("Payoff_A = P(CC)*V(A,CC) + P(CD)*V(A,CD) + P(DC)*V(A,DC) + P(DD)*V(A,DD)")
    print(f"Payoff_A = ({probs_QQ[0]:.2f} * {payoff_values['A'][0]}) + ({probs_QQ[1]:.2f} * {payoff_values['A'][1]}) + ({probs_QQ[2]:.2f} * {payoff_values['A'][2]}) + ({probs_QQ[3]:.2f} * {payoff_values['A'][3]}) = {payoff_A_QQ:.2f}")

    print("\nBob's Equilibrium Payoff Calculation:")
    print("Payoff_B = P(CC)*V(B,CC) + P(CD)*V(B,CD) + P(DC)*V(B,DC) + P(DD)*V(B,DD)")
    print(f"Payoff_B = ({probs_QQ[0]:.2f} * {payoff_values['B'][0]}) + ({probs_QQ[1]:.2f} * {payoff_values['B'][1]}) + ({probs_QQ[2]:.2f} * {payoff_values['B'][2]}) + ({probs_QQ[3]:.2f} * {payoff_values['B'][3]}) = {payoff_B_QQ:.2f}")

solve_quantum_prisoners_dilemma()
<<<5>>>