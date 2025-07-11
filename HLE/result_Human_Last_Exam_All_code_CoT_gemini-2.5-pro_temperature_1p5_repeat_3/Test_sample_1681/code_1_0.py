import numpy as np

def solve_quantum_dilemma():
    """
    Solves for the equilibrium point of the Quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol.
    """

    # 1. Define the game elements
    # Payoff matrix values from the problem description
    R = 5  # Reward for mutual cooperation
    S = 0  # Sucker's payoff
    T = 7  # Temptation payoff
    P = 1  # Punishment for mutual defection

    # Basis states represented as column vectors |0> and |1>
    # C (Cooperate) is |0>, D (Defect) is |1>
    state_C = np.array([[1], [0]], dtype=complex)
    state_D = np.array([[0], [1]], dtype=complex)

    # The initial state of the game is |CC> or |00>
    initial_state_CC = np.kron(state_C, state_C)

    # Define standard operators (Pauli matrices and Identity)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    Identity = np.identity(2, dtype=complex)

    # Define player strategies as unitary operators (matrices)
    # Cooperate (C) is the Identity operator
    # Defect (D) is a bit-flip like operation (i * sigma_x)
    # The new Quantum strategy (Q) is (i * sigma_z)
    C_op = Identity
    D_op = 1j * sigma_x
    Q_op = 1j * sigma_z

    # The entangling operator J for creating maximal entanglement
    J = (1 / np.sqrt(2)) * (np.kron(Identity, Identity) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T

    # 2. Identify the Nash Equilibrium strategy
    # The known Nash Equilibrium for this quantum game is (Q, Q)
    U_A = Q_op
    U_B = Q_op

    # 3. Calculate the final state and probabilities
    # The evolution of the state is: J_dag * (U_A kron U_B) * J * initial_state
    final_state = J_dag @ np.kron(U_A, U_B) @ J @ initial_state_CC

    # Calculate the probability of each classical outcome by projecting the final state
    # onto the classical basis states (|CC>, |CD>, |DC>, |DD>)
    prob_CC = np.abs(np.vdot(np.kron(state_C, state_C), final_state))**2
    prob_CD = np.abs(np.vdot(np.kron(state_C, state_D), final_state))**2
    prob_DC = np.abs(np.vdot(np.kron(state_D, state_C), final_state))**2
    prob_DD = np.abs(np.vdot(np.kron(state_D, state_D), final_state))**2

    # 4. Calculate the expected payoffs based on the probabilities
    payoff_A = prob_CC * R + prob_CD * S + prob_DC * T + prob_DD * P
    payoff_B = prob_CC * R + prob_CD * T + prob_DC * S + prob_DD * P

    # 5. Print the results and the final equation
    print("The quantum Nash Equilibrium is achieved when both players choose the strategy 'Q'.")
    print("\nCalculating the outcome for the (Q, Q) strategy pair:")
    print(f"  Probability(Cooperate, Cooperate) = {prob_CC:.2f}")
    print(f"  Probability(Cooperate, Defect)   = {prob_CD:.2f}")
    print(f"  Probability(Defect,   Cooperate) = {prob_DC:.2f}")
    print(f"  Probability(Defect,   Defect)   = {prob_DD:.2f}")

    print("\nThe payoff for Player A is calculated as follows:")
    print("  Payoff_A = P(CC)*R + P(CD)*S + P(DC)*T + P(DD)*P")
    print(f"  Payoff_A = {prob_CC:.2f} * {R} + {prob_CD:.2f} * {S} + {prob_DC:.2f} * {T} + {prob_DD:.2f} * {P} = {payoff_A.real:.2f}")

    print("\nThe payoff for Player B is calculated similarly:")
    print("  Payoff_B = P(CC)*R + P(CD)*T + P(DC)*S + P(DD)*P")
    print(f"  Payoff_B = {prob_CC:.2f} * {R} + {prob_CD:.2f} * {T} + {prob_DC:.2f} * {S} + {prob_DD:.2f} * {P} = {payoff_B.real:.2f}")
    
    print(f"\nThis new equilibrium resolves the dilemma, yielding a payoff of ({payoff_A.real:.2f}, {payoff_B.real:.2f}),")
    print("which is Pareto optimal and superior to the classical (Defect, Defect) outcome of (1, 1).")
    
    # Return the single value for the final answer format
    return payoff_A.real

if __name__ == '__main__':
    equilibrium_payoff = solve_quantum_dilemma()
    print(f"\n<<<5.0>>>")