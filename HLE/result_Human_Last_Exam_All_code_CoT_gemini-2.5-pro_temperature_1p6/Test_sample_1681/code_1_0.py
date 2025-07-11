import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates and explains the equilibrium point for the quantum Prisoner's Dilemma.
    """
    # 1. SETUP THE QUANTUM GAME

    # Payoff matrix from the problem description
    # P[i][j] = (Payoff_A, Payoff_B) for outcome row i, col j
    # (C,C): (5,5), (C,D): (0,7), (D,C): (7,0), (D,D): (1,1)
    payoffs_A = np.array([[5, 0], [7, 1]])
    payoffs_B = np.array([[5, 7], [0, 1]])

    # Initial state |00> (C,C)
    psi_0 = np.array([1, 0, 0, 0], dtype=complex)

    # Entangling Operator J and its conjugate transpose J_dag
    # J = exp(i * gamma/2 * sigma_x ⊗ sigma_x), with gamma=pi/2 for max entanglement
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T

    # Player Strategies as Unitary Operators from SU(2)
    # C = U(0,0), D = U(pi,0), Q = U(pi/2, pi/2)
    # Classical Cooperate (Identity)
    C = np.array([[1, 0], [0, 1]], dtype=complex)
    # Classical Defect
    D = np.array([[0, 1], [-1, 0]], dtype=complex)
    # Eisert's "Miracle Move" Quantum Strategy
    Q = np.array([[1j/np.sqrt(2), 1/np.sqrt(2)], [-1/np.sqrt(2), -1j/np.sqrt(2)]], dtype=complex)

    def calculate_payoffs(U_A, U_B):
        """Calculates final state, probabilities, and payoffs for strategies U_A, U_B."""
        # The full operator applied to the initial state
        game_operator = J_dag @ np.kron(U_A, U_B) @ J
        
        # Calculate the final state vector
        psi_final = game_operator @ psi_0
        
        # Probabilities of the 4 classical outcomes (CC, CD, DC, DD)
        probabilities = np.abs(psi_final)**2
        
        # Expected payoffs
        payoff_A = np.sum(probabilities * payoffs_A.flatten())
        payoff_B = np.sum(probabilities * payoffs_B.flatten())
        
        return payoff_A, payoff_B, probabilities

    # 2. IDENTIFY AND VERIFY THE EQUILIBRIUM
    print("--- Analysis of Quantum Prisoner's Dilemma Equilibrium ---")
    print("\nPayoff Matrix (Player A, Player B):")
    print(f"  Cooperate(B) Defect(B)")
    print(f"C(A)   ({payoffs_A[0,0]},{payoffs_B[0,0]})     ({payoffs_A[0,1]},{payoffs_B[0,1]})")
    print(f"D(A)   ({payoffs_A[1,0]},{payoffs_B[1,0]})      ({payoffs_A[1,1]},{payoffs_B[1,1]})")
    
    # --- Case 1: The (Q, Q) Equilibrium ---
    print("\n1. Calculating the Equilibrium Point with Strategy (Q, Q):")
    pA_QQ, pB_QQ, probs_QQ = calculate_payoffs(Q, Q)
    
    print(f"\nFinal probabilities for outcomes (CC, CD, DC, DD) are: [{probs_QQ[0]:.2f}, {probs_QQ[1]:.2f}, {probs_QQ[2]:.2f}, {probs_QQ[3]:.2f}]")
    
    print("\nPayoff Calculation:")
    print(f"Player A = ({probs_QQ[0]:.2f} * {payoffs_A[0,0]}) + ({probs_QQ[1]:.2f} * {payoffs_A[0,1]}) + ({probs_QQ[2]:.2f} * {payoffs_A[1,0]}) + ({probs_QQ[3]:.2f} * {payoffs_A[1,1]}) = {pA_QQ:.2f}")
    print(f"Player B = ({probs_QQ[0]:.2f} * {payoffs_B[0,0]}) + ({probs_QQ[1]:.2f} * {payoffs_B[0,1]}) + ({probs_QQ[2]:.2f} * {payoffs_B[1,0]}) + ({probs_QQ[3]:.2f} * {payoffs_B[1,1]}) = {pB_QQ:.2f}")
    print(f"\nThe equilibrium payoff is ({pA_QQ:.2f}, {pB_QQ:.2f}), avoiding the classical dilemma.")

    # --- Case 2: Verify it's a Nash Equilibrium ---
    print("\n\n2. Verifying if (Q, Q) is a Nash Equilibrium:")
    print("We check if Player B can get a better payoff by deviating while Player A plays Q.")

    # Player B deviates to Cooperate (C)
    pA_QC, pB_QC, _ = calculate_payoffs(Q, C)
    print(f"\nIf Player B deviates to 'C', the payoffs (A, B) become: ({pA_QC:.2f}, {pB_QC:.2f})")

    # Player B deviates to Defect (D)
    pA_QD, pB_QD, _ = calculate_payoffs(Q, D)
    print(f"If Player B deviates to 'D', the payoffs (A, B) become: ({pA_QD:.2f}, {pB_QD:.2f})")
    
    print("\nComparison of Player B's payoffs:")
    print(f"  - Payoff from staying with Q: {pB_QQ:.2f}")
    print(f"  - Payoff from deviating to C: {pB_QC:.2f}")
    print(f"  - Payoff from deviating to D: {pB_QD:.2f}")

    if pB_QQ >= pB_QC and pB_QQ >= pB_QD:
        print("\nConclusion: Player B has no incentive to unilaterally deviate from Q.")
        print("The strategy profile (Q, Q) is a stable Nash Equilibrium.")
    else:
        print("\nConclusion: (Q, Q) is not a Nash Equilibrium.")


solve_quantum_prisoners_dilemma()