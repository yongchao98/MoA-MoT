import numpy as np

def run_quantum_pd_analysis():
    """
    Analyzes the Prisoner's Dilemma in a quantum setting to find the equilibrium point.
    """
    # Define Pauli matrices and Identity
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    # Payoff values for Alice and Bob for outcomes |00>, |01>, |10>, |11>
    # Outcomes correspond to (C,C), (C,D), (D,C), (D,D)
    # Alice's payoffs: Cooperate (C) is row 1, Defect (D) is row 2
    # Bob's payoffs: Cooperate (C) is col 1, Defect (D) is col 2
    # Matrix: [[(5,5), (0,7)], [(7,0), (1,1)]]
    payoffs_A = np.array([5, 0, 7, 1])
    payoffs_B = np.array([5, 7, 0, 1])
    
    # Define strategy operators
    C = I      # Cooperate
    D = sx     # Defect
    Q = -1j * sz # Quantum strategy

    # Define the entangling operator J
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sx, sx))
    J_dag = J.conj().T

    # Initial state |00>
    psi_0 = np.array([1, 0, 0, 0]).reshape(4, 1)

    def calculate_payoffs(U_A, U_B, name_A, name_B):
        """Calculates and prints the outcome of a game with given strategies."""
        print(f"\n--- Calculating Payoff for strategy ({name_A}, {name_B}) ---")
        
        # Combine strategies using Kronecker product
        U_total = np.kron(U_A, U_B)
        
        # Calculate the final state evolution operator and the final state
        operator = J_dag @ U_total @ J
        psi_f = operator @ psi_0
        
        # Calculate probabilities of measuring each classical outcome
        probs = (np.abs(psi_f)**2).flatten()

        # Calculate expected payoffs
        payoff_A = np.sum(probs * payoffs_A)
        payoff_B = np.sum(probs * payoffs_B)
        
        print(f"Probabilities [CC, CD, DC, DD]: "
              f"[{probs[0]:.2f}, {probs[1]:.2f}, {probs[2]:.2f}, {probs[3]:.2f}]")
        
        print("Final payoff equation for Alice:")
        print(f"Payoff_A = ({payoffs_A[0]})*({probs[0]:.2f}) + ({payoffs_A[1]})*({probs[1]:.2f}) + "
              f"({payoffs_A[2]})*({probs[2]:.2f}) + ({payoffs_A[3]})*({probs[3]:.2f}) = {payoff_A:.2f}")

        print("Final payoff equation for Bob:")
        print(f"Payoff_B = ({payoffs_B[0]})*({probs[0]:.2f}) + ({payoffs_B[1]})*({probs[1]:.2f}) + "
              f"({payoffs_B[2]})*({probs[2]:.2f}) + ({payoffs_B[3]})*({probs[3]:.2f}) = {payoff_B:.2f}")

        print(f"Final Payoffs (Alice, Bob): ({payoff_A:.2f}, {payoff_B:.2f})")
        return payoff_A, payoff_B

    print("Analyzing the quantum Prisoner's Dilemma equilibrium.")
    print("Strategies: C (Cooperate), D (Defect), Q (Quantum).")

    # 1. Test the (Q, Q) strategy profile - the potential new equilibrium
    q_q_payoffs = calculate_payoffs(Q, Q, "Q", "Q")

    print("\nNow, let's check if Alice has an incentive to deviate while Bob plays Q.")
    
    # 2. Alice deviates to Defect (D)
    d_q_payoffs = calculate_payoffs(D, Q, "D", "Q")
    print(f"Result: Alice's payoff by deviating to D ({d_q_payoffs[0]:.2f}) is less than her payoff from (Q,Q) ({q_q_payoffs[0]:.2f}).")
    
    # 3. Alice deviates to Cooperate (C)
    c_q_payoffs = calculate_payoffs(C, Q, "C", "Q")
    print(f"Result: Alice's payoff by deviating to C ({c_q_payoffs[0]:.2f}) is less than her payoff from (Q,Q) ({q_q_payoffs[0]:.2f}).")

    print("\n--- Conclusion ---")
    print("Since neither player can improve their payoff by unilaterally changing their strategy from Q,")
    print("(Q, Q) is a Nash Equilibrium. This equilibrium yields the payoff (5.00, 5.00), which corresponds")
    print("to the optimal outcome of mutual cooperation, thus resolving the dilemma.")
    print("The equilibrium point is the payoff profile (5, 5).")

run_quantum_pd_analysis()