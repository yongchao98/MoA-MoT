import numpy as np

def quantum_prisoners_dilemma_equilibrium():
    """
    Calculates and prints the payoff for the quantum equilibrium point
    in the Prisoner's Dilemma using the EWL protocol.
    """
    # 1. Payoff matrix values
    # (Cooperate, Cooperate), (Cooperate, Defect), (Defect, Cooperate), (Defect, Defect)
    R, S, T, P = 5, 0, 7, 1
    
    # 2. Define fundamental quantum objects
    # Identity and Pauli matrices
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    # Basis vectors |0> (Cooperate) and |1> (Defect)
    ket0 = np.array([[1], [0]], dtype=complex)
    ket1 = np.array([[0], [1]], dtype=complex)

    # 3. Define the players' strategies as unitary operators
    # Classical Cooperate: Identity operator
    C = I
    # Quantum "Miracle" Strategy Q
    Q = -1j * sz

    # 4. Set up the game using the EWL protocol
    # Entangling operator J for maximal entanglement
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sx, sx))
    # Disentangling operator is the conjugate transpose of J
    J_dag = J.conj().T

    # Initial state is |CC> or |00>
    psi_0 = np.kron(ket0, ket0)
    
    # The arbiter entangles the initial state
    psi_initial = J @ psi_0

    # 5. Calculate the outcome for the (Q, Q) equilibrium
    # Both players apply the quantum strategy Q
    U_A = Q
    U_B = Q
    
    # The final state is calculated by applying player strategies and then disentangling
    psi_final = J_dag @ np.kron(U_A, U_B) @ psi_initial

    # 6. Calculate outcome probabilities by projecting the final state
    # onto the classical basis states
    p_cc = np.abs((np.kron(ket0, ket0).conj().T @ psi_final)[0, 0])**2
    p_cd = np.abs((np.kron(ket0, ket1).conj().T @ psi_final)[0, 0])**2
    p_dc = np.abs((np.kron(ket1, ket0).conj().T @ psi_final)[0, 0])**2
    p_dd = np.abs((np.kron(ket1, ket1).conj().T @ psi_final)[0, 0])**2

    # 7. Calculate the payoffs for each player
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd
    
    # 8. Print the result
    # We round the results to remove potential floating point inaccuracies
    r, s, t, p = R, S, T, P
    pcc, pcd, pdc, pdd = round(p_cc), round(p_cd), round(p_dc), round(p_dd)
    pa, pb = round(payoff_A), round(payoff_B)

    print("Quantum Prisoner's Dilemma Equilibrium (Q,Q):")
    print(f"Final state probabilities: P(CC)={pcc}, P(CD)={pcd}, P(DC)={pdc}, P(DD)={pdd}")
    print("\nFinal Payoff Equation for Player A:")
    print(f"Payoff(A) = P(CC)*{r} + P(CD)*{s} + P(DC)*{t} + P(DD)*{p}")
    print(f"Payoff(A) = {pcc}*{r} + {pcd}*{s} + {pdc}*{t} + {pdd}*{p} = {pa}")

    print("\nFinal Payoff Equation for Player B:")
    print(f"Payoff(B) = P(CC)*{r} + P(CD)*{t} + P(DC)*{s} + P(DD)*{p}")
    print(f"Payoff(B) = {pcc}*{r} + {pcd}*{t} + {pdc}*{s} + {pdd}*{p} = {pb}")

    print(f"\nThe equilibrium point is the strategy pair (Q, Q), with a payoff of ({pa}, {pb}).")

if __name__ == '__main__':
    quantum_prisoners_dilemma_equilibrium()