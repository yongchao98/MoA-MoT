import numpy as np

def run_quantum_prisoners_dilemma_equilibrium():
    """
    Calculates and explains the equilibrium point for the quantum Prisoner's Dilemma.

    In the quantum version, players' states are entangled, and their strategies are
    unitary operators. A new strategy, 'Q', becomes available. A new Nash Equilibrium
    (Q, Q) emerges, leading to a Pareto-optimal outcome.
    """

    # --- 1. Define Quantum Primitives ---
    # Basis states |0> (Cooperate) and |1> (Defect)
    q0 = np.array([[1], [0]], dtype=complex)  # |0> (Cooperate)
    q1 = np.array([[0], [1]], dtype=complex)  # |1> (Defect)

    # Basis for the 2-qubit system
    psi00 = np.kron(q0, q0)  # |00> or |CC>

    # Pauli matrices and Identity
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # --- 2. Define Game Setup ---
    # Payoff values [PlayerA, PlayerB]
    # (Cooperate, Cooperate), (Cooperate, Defect), (Defect, Cooperate), (Defect, Defect)
    payoffs = {
        "CC": (5, 5), "CD": (0, 7), "DC": (7, 0), "DD": (1, 1)
    }

    # Entangling operator J for maximal entanglement (gamma = pi/2)
    gamma = np.pi / 2
    J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(sigma_x, sigma_x)
    J_dag = J.conj().T  # Adjoint (disentangling operator)

    # Initial entangled state
    initial_state = J @ psi00

    # --- 3. Define Player Strategies ---
    # C = Cooperate, D = Defect, Q = Quantum
    U_C = I
    U_D = sigma_x
    U_Q = 1j * sigma_z # The "miracle" quantum move

    # --- 4. Calculate the Outcome for the Equilibrium Strategy (Q, Q) ---
    # At equilibrium, both players choose the quantum strategy Q.
    U_A = U_Q
    U_B = U_Q

    # The combined strategy operator
    U_total = np.kron(U_A, U_B)

    # The final state of the game after strategies and disentanglement
    final_state = J_dag @ U_total @ initial_state

    # Probabilities are the squared magnitudes of the final state vector's components
    # The components correspond to |00>, |01>, |10>, |11>
    prob_00 = np.abs(final_state[0, 0])**2 # P(Cooperate, Cooperate)
    prob_01 = np.abs(final_state[1, 0])**2 # P(Cooperate, Defect)
    prob_10 = np.abs(final_state[2, 0])**2 # P(Defect, Cooperate)
    prob_11 = np.abs(final_state[3, 0])**2 # P(Defect, Defect)

    # --- 5. Calculate and Print the Final Payoffs ---
    # Expected payoff is the sum of (payoff * probability) for each outcome
    payoff_A = (payoffs["CC"][0] * prob_00 +
                payoffs["CD"][0] * prob_01 +
                payoffs["DC"][0] * prob_10 +
                payoffs["DD"][0] * prob_11)

    payoff_B = (payoffs["CC"][1] * prob_00 +
                payoffs["CD"][1] * prob_01 +
                payoffs["DC"][1] * prob_10 +
                payoffs["DD"][1] * prob_11)

    print("The quantum equilibrium strategy is (Q, Q), where Q is a special quantum move.")
    print("This leads to the following outcome probabilities:")
    print(f"P(Cooperate, Cooperate) = {prob_00:.2f}")
    print(f"P(Cooperate, Defect)   = {prob_01:.2f}")
    print(f"P(Defect,   Cooperate) = {prob_10:.2f}")
    print(f"P(Defect,   Defect)   = {prob_11:.2f}\n")

    print("The equilibrium payoff for each player is calculated as:")
    print("Payoff = P(CC)*V(CC) + P(CD)*V(CD) + P(DC)*V(DC) + P(DD)*V(DD)\n")

    print("Player A's payoff calculation:")
    print(f"{payoff_A:.2f} = {payoffs['CC'][0]} * {prob_00:.2f} + {payoffs['CD'][0]} * {prob_01:.2f} + {payoffs['DC'][0]} * {prob_10:.2f} + {payoffs['DD'][0]} * {prob_11:.2f}")

    print("\nPlayer B's payoff calculation:")
    print(f"{payoff_B:.2f} = {payoffs['CC'][1]} * {prob_00:.2f} + {payoffs['CD'][1]} * {prob_01:.2f} + {payoffs['DC'][1]} * {prob_10:.2f} + {payoffs['DD'][1]} * {prob_11:.2f}")

    print(f"\nThus, the equilibrium point is a payoff of ({int(payoff_A)}, {int(payoff_B)}), resolving the dilemma.")
    
# Execute the function
run_quantum_prisoners_dilemma_equilibrium()