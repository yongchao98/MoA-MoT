import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point payoff for the Quantum Prisoner's Dilemma.
    """
    # 1. Define the payoff matrix values
    # (Cooperate, Cooperate): (R, R)
    # (Cooperate, Defect):   (S, T)
    # (Defect,   Cooperate): (T, S)
    # (Defect,   Defect):   (P, P)
    R = 5  # Reward
    S = 0  # Sucker
    T = 7  # Temptation
    P = 1  # Punishment

    # 2. Define quantum operators and basis states
    # Basis states |C> (Cooperate) and |D> (Defect)
    C = np.array([[1], [0]])
    D = np.array([[0], [1]])

    # Identity and Pauli Matrices
    I = np.eye(2, dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # The Quantum strategy 'Q'
    Q = 1j * sigma_z

    # 3. Set up the game according to the EWL protocol
    # The entangling operator J for maximal entanglement (gamma = pi/2)
    # J = exp(i * gamma/2 * (sigma_x cross sigma_x))
    # Using cos(pi/4)*I⊗I + i*sin(pi/4)*σx⊗σx gives the same J as below
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(np.array([[0,1],[1,0]]), np.array([[0,1],[1,0]])))
    J_dag = J.conj().T

    # The initial state is |CC> = |C> ⊗ |C>
    psi_0 = np.kron(C, C)

    # 4. Calculate the final state for the equilibrium strategy (Q, Q)
    # The operator for the players' moves
    U_game = np.kron(Q, Q)
    
    # The evolution of the state vector
    psi_final = J_dag @ U_game @ J @ psi_0

    # 5. Calculate outcome probabilities
    # Define the four classical basis vectors for the 2-qubit system
    basis_CC = np.kron(C, C) # Alice Cooperates, Bob Cooperates
    basis_CD = np.kron(C, D) # Alice Cooperates, Bob Defects
    basis_DC = np.kron(D, C) # Alice Defects,   Bob Cooperates
    basis_DD = np.kron(D, D) # Alice Defects,   Bob Defects

    # Probabilities are the squared magnitudes of the final state's components
    prob_CC = np.abs(basis_CC.conj().T @ psi_final)**2
    prob_CD = np.abs(basis_CD.conj().T @ psi_final)**2
    prob_DC = np.abs(basis_DC.conj().T @ psi_final)**2
    prob_DD = np.abs(basis_DD.conj().T @ psi_final)**2
    
    # 6. Calculate the expected payoffs
    payoff_A = prob_CC * R + prob_CD * S + prob_DC * T + prob_DD * P
    payoff_B = prob_CC * R + prob_CD * T + prob_DC * S + prob_DD * P

    # 7. Print the results
    # We round the results to handle potential floating point inaccuracies.
    prob_CC = np.round(prob_CC.item(), 4)
    prob_CD = np.round(prob_CD.item(), 4)
    prob_DC = np.round(prob_DC.item(), 4)
    prob_DD = np.round(prob_DD.item(), 4)
    payoff_A = np.round(payoff_A.item(), 4)
    payoff_B = np.round(payoff_B.item(), 4)

    print("Quantum Prisoner's Dilemma Equilibrium Analysis:")
    print("-" * 50)
    print("The equilibrium strategy profile is (Q, Q).")
    print(f"Outcome Probabilities for (Q, Q):")
    print(f"  P(Cooperate, Cooperate) = {prob_CC}")
    print(f"  P(Cooperate, Defect)   = {prob_CD}")
    print(f"  P(Defect,   Cooperate) = {prob_DC}")
    print(f"  P(Defect,   Defect)   = {prob_DD}")
    print("-" * 50)
    print("Expected Payoff Calculation for Player A:")
    print(f"Payoff(A) = P(CC)*R + P(CD)*S + P(DC)*T + P(DD)*P")
    print(f"Payoff(A) = {prob_CC}*{R} + {prob_CD}*{S} + {prob_DC}*{T} + {prob_DD}*{P} = {payoff_A}")
    print("\nExpected Payoff Calculation for Player B:")
    print(f"Payoff(B) = P(CC)*R + P(CD)*T + P(DC)*S + P(DD)*P")
    print(f"Payoff(B) = {prob_CC}*{R} + {prob_CD}*{T} + {prob_DC}*{S} + {prob_DD}*{P} = {payoff_B}")
    print("-" * 50)
    print(f"The equilibrium point results in the payoff: ({payoff_A}, {payoff_B})")

    # Final answer format for one player's payoff
    return payoff_A

# Execute the function and capture the final answer
equilibrium_payoff = solve_quantum_prisoners_dilemma()
print(f"\n<<<({equilibrium_payoff}, {equilibrium_payoff})>>>")
