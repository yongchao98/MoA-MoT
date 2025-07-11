import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol with maximal entanglement.
    """
    # 1. Define the payoff matrix from the problem description.
    # Payoffs for (CC, CD, DC, DD) outcomes.
    # R=Reward, S=Sucker, T=Temptation, P=Punishment
    R, S, T, P = 5, 0, 7, 1
    payoffs_A = np.array([R, S, T, P])
    payoffs_B = np.array([R, T, S, P])

    # 2. Set up the quantum game framework (EWL Protocol).
    # We use a maximally entangling operator J, corresponding to gamma = pi/2.
    # J = (1/sqrt(2)) * (I ⊗ I + i * σ_x ⊗ σ_x)
    I4 = np.identity(4, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_x_kron_sigma_x = np.kron(sigma_x, sigma_x)
    
    J = (1 / np.sqrt(2)) * (I4 + 1j * sigma_x_kron_sigma_x)
    J_dag = J.conj().T  # Hermitian conjugate of J

    # 3. The initial state is |CC>, which is |00> in quantum computing.
    # psi_0 = [1, 0, 0, 0]
    psi_0 = np.zeros(4, dtype=complex)
    psi_0[0] = 1

    # 4. Define the quantum strategy 'Q' that forms the new equilibrium.
    # This strategy is a unitary operator: Q = [[i, 0], [0, -i]]
    Q = np.array([[1j, 0], [0, -1j]])

    # Both players choose the strategy Q.
    U_A = Q
    U_B = Q
    U_kron_U = np.kron(U_A, U_B)

    # 5. Simulate the game's evolution to find the final state.
    # final_state = J_dag * (U_A ⊗ U_B) * J * initial_state
    final_state = J_dag @ U_kron_U @ J @ psi_0

    # 6. Calculate outcome probabilities from the final state vector.
    # Probabilities are the squared magnitudes of the final state's amplitudes.
    probs = np.abs(final_state)**2

    # 7. Calculate the expected payoff for each player.
    payoff_A = np.sum(probs * payoffs_A)
    payoff_B = np.sum(probs * payoffs_B)

    # 8. Print the results and the calculation.
    print("In the Quantum Prisoner's Dilemma with maximal entanglement, a new Nash Equilibrium emerges.")
    print("This equilibrium corresponds to both players choosing the quantum strategy Q.\n")
    print("The final probabilities of the classical outcomes (CC, CD, DC, DD) are:")
    print(f"P(CC)={probs[0]:.2f}, P(CD)={probs[1]:.2f}, P(DC)={probs[2]:.2f}, P(DD)={probs[3]:.2f}\n")

    print("Player A's Payoff Calculation:")
    print(f"E_A = P(CC)xR + P(CD)xS + P(DC)xT + P(DD)xP")
    print(f"E_A = {probs[0]:.2f} x {R} + {probs[1]:.2f} x {S} + {probs[2]:.2f} x {T} + {probs[3]:.2f} x {P} = {payoff_A:.0f}\n")
    
    print("Player B's Payoff Calculation:")
    print(f"E_B = P(CC)xR + P(CD)xT + P(DC)xS + P(DD)xP")
    print(f"E_B = {probs[0]:.2f} x {R} + {probs[1]:.2f} x {T} + {probs[2]:.2f} x {S} + {probs[3]:.2f} x {P} = {payoff_B:.0f}\n")
    
    print(f"The equilibrium point (Payoff A, Payoff B) is: ({payoff_A:.0f}, {payoff_B:.0f})")

solve_quantum_prisoners_dilemma()