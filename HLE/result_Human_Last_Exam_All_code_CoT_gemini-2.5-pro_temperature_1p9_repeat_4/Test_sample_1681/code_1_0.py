import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the quantum Prisoner's Dilemma.

    This function implements the Eisert-Wilkens-Lewenstein (EWL) protocol
    with maximal entanglement to find the payoff at the quantum Nash Equilibrium.
    """

    # Payoff matrix values
    # (Cooperate, Cooperate): R=5
    # (Cooperate, Defect): S=0 (for Player 1), T=7 (for Player 2)
    # (Defect, Cooperate): T=7 (for Player 1), S=0 (for Player 2)
    # (Defect, Defect): P=1
    R, S, T, P = 5, 0, 7, 1

    # Define basis states |00>, |01>, |10>, |11>
    # |0> = Cooperate, |1> = Defect
    psi_00 = np.array([1, 0, 0, 0])  # Represents initial state |CC>

    # Pauli-X and Identity matrices for constructing the entangler
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])

    # Entangling operator J for maximal entanglement (gamma = pi/2)
    # J = exp(i * pi/4 * sigma_x @ sigma_x) = (I @ I + i * sigma_x @ sigma_x) / sqrt(2)
    # In numpy, '@' is matrix multiplication, `np.kron` is the tensor product.
    J = (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x)) / np.sqrt(2)
    J_dag = J.conj().T

    # The Quantum Nash Equilibrium strategy 'Q'
    # This is a known strategy that forms a new equilibrium in the quantum game.
    Q = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)

    # Both players apply the quantum strategy Q
    Q_otimes_Q = np.kron(Q, Q)

    # Calculate the final state before measurement
    # psi_final = J_dag * (Q @ Q) * J * psi_00
    psi_mid = J @ psi_00
    psi_transformed = Q_otimes_Q @ psi_mid
    psi_final = J_dag @ psi_transformed
    
    # Calculate probabilities of the four classical outcomes
    p_cc = np.abs(psi_final[0])**2
    p_cd = np.abs(psi_final[1])**2
    p_dc = np.abs(psi_final[2])**2
    p_dd = np.abs(psi_final[3])**2

    # Calculate the expected payoffs for each player
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd
    
    # Print the equation and the final payoff for Player A
    # The game is symmetric, so Player B's payoff is the same.
    print("At the quantum Nash Equilibrium (Q, Q):")
    print(f"Probabilities: P(CC)={p_cc:.2f}, P(CD)={p_cd:.2f}, P(DC)={p_dc:.2f}, P(DD)={p_dd:.2f}")
    print("\nFinal Payoff Calculation:")
    print(f"Payoff(A) = ({R} * {p_cc:.2f}) + ({S} * {p_cd:.2f}) + ({T} * {p_dc:.2f}) + ({P} * {p_dd:.2f}) = {payoff_A:.2f}")
    print(f"Payoff(B) = ({R} * {p_cc:.2f}) + ({T} * {p_cd:.2f}) + ({S} * {p_dc:.2f}) + ({P} * {p_dd:.2f}) = {payoff_B:.2f}")
    
    # Return the equilibrium point value
    # For a symmetric game, this is the payoff each player receives.
    return payoff_A

# Run the calculation and store the result for the final answer
equilibrium_payoff = solve_quantum_prisoners_dilemma()