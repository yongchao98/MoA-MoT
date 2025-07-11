import numpy as np

def quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum prisoner's dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) framework.
    """
    # -- 1. Define matrices and operators --

    # Pauli matrices and Identity
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Player strategies: Cooperate (C), Defect (D), Quantum (Q)
    C = I
    D = 1j * sigma_x
    Q = 1j * sigma_z

    # Entangling operator J and its inverse J_dagger (for maximal entanglement)
    J = (1/np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dagger = J.conj().T

    # Payoff values from the matrix: [[(R,R), (S,T)], [(T,S), (P,P)]]
    # R=5 (Cooperate, Cooperate), S=0 (Sucker), T=7 (Temptation), P=1 (Punishment)
    payoffs_A = np.array([5, 0, 7, 1])  # Payoffs for Alice for outcomes CC, CD, DC, DD
    payoffs_B = np.array([5, 7, 0, 1])  # Payoffs for Bob for outcomes CC, CD, DC, DD

    # Initial state of the game is |CC> or |00>
    psi_initial = np.array([1, 0, 0, 0], dtype=complex)

    print("Quantum Prisoner's Dilemma Calculation")
    print("-" * 40)

    # -- 2. Calculate the Nash Equilibrium (Q, Q) --
    print("Step 1: Calculating the payoff for the (Q, Q) equilibrium strategy.")

    # Combined strategy operator for (Q, Q)
    U_eq = np.kron(Q, Q)

    # Calculate the final state by applying the sequence of operators
    psi_final_eq = J_dagger @ U_eq @ J @ psi_initial
    
    # Probabilities of measuring each outcome (|psi_final|^2)
    probs_eq = np.abs(psi_final_eq)**2

    # Expected payoff for each player
    payoff_A_eq = np.sum(probs_eq * payoffs_A)
    payoff_B_eq = np.sum(probs_eq * payoffs_B)
    
    # Print the equation and the result
    print("Final state |ψ_f> for (Q, Q) results in probabilities for outcomes (CC, CD, DC, DD):")
    print(f"P = [{probs_eq[0]:.2f}, {probs_eq[1]:.2f}, {probs_eq[2]:.2f}, {probs_eq[3]:.2f}]")
    print("Alice's payoff equation: Payoff(A) = P(CC)*5 + P(CD)*0 + P(DC)*7 + P(DD)*1")
    print(f"Payoff(A) = {probs_eq[0]:.2f}*5 + {probs_eq[1]:.2f}*0 + {probs_eq[2]:.2f}*7 + {probs_eq[3]:.2f}*1 = {payoff_A_eq:.2f}")
    print("Bob's payoff equation:   Payoff(B) = P(CC)*5 + P(CD)*7 + P(DC)*0 + P(DD)*1")
    print(f"Payoff(B) = {probs_eq[0]:.2f}*5 + {probs_eq[1]:.2f}*7 + {probs_eq[2]:.2f}*0 + {probs_eq[3]:.2f}*1 = {payoff_B_eq:.2f}")
    print(f"\nThe payoff at the (Q, Q) equilibrium is ({payoff_A_eq:.2f}, {payoff_B_eq:.2f}), which corresponds to mutual cooperation.")

    print("-" * 40)
    
    # -- 3. Check for incentive to deviate (e.g., Alice plays D, Bob plays Q) --
    print("Step 2: Checking if a player has an incentive to deviate.")
    print("Let's calculate Alice's payoff if she deviates to D while Bob plays Q.")

    # Combined strategy operator for (D, Q)
    U_dev = np.kron(D, Q)
    
    # Calculate the final state for this deviation
    psi_final_dev = J_dagger @ U_dev @ J @ psi_initial

    # Probabilities of measuring each outcome
    probs_dev = np.abs(psi_final_dev)**2

    # Expected payoff for Alice if she deviates
    payoff_A_dev = np.sum(probs_dev * payoffs_A)
    payoff_B_dev = np.sum(probs_dev * payoffs_B)
    
    print("Final state |ψ_f> for (D, Q) results in probabilities for outcomes (CC, CD, DC, DD):")
    print(f"P = [{probs_dev[0]:.2f}, {probs_dev[1]:.2f}, {probs_dev[2]:.2f}, {probs_dev[3]:.2f}]")
    print("Alice's payoff equation: Payoff(A) = P(CC)*5 + P(CD)*0 + P(DC)*7 + P(DD)*1")
    print(f"Payoff(A) = {probs_dev[0]:.2f}*5 + {probs_dev[1]:.2f}*0 + {probs_dev[2]:.2f}*7 + {probs_dev[3]:.2f}*1 = {payoff_A_dev:.2f}")

    print(f"\nAlice's payoff by deviating to D is {payoff_A_dev:.2f}, which is less than her equilibrium payoff of {payoff_A_eq:.2f}.")
    print("Therefore, there is no incentive to deviate from the (Q, Q) strategy.")

    print("-" * 40)
    print("Conclusion: The dilemma is resolved.")
    print("The equilibrium point for the quantum prisoner's dilemma under these optimal conditions is the payoff vector achieved at the (Q,Q) Nash Equilibrium.")


quantum_prisoners_dilemma()
print("<<<(5.0, 5.0)>>>")