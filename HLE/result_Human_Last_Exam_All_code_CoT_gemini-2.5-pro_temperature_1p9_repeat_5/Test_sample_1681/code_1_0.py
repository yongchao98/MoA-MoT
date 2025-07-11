import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates and prints the payoff for the equilibrium point in a
    quantum Prisoner's Dilemma with a specific payoff matrix.
    """
    # Payoff matrix values from the problem description
    # (Cooperate, Cooperate): (R, R) -> (5, 5)
    # (Cooperate, Defect):   (S, T) -> (0, 7)
    # (Defect, Cooperate):   (T, S) -> (7, 0)
    # (Defect, Defect):      (P, P) -> (1, 1)
    R = 5  # Reward for mutual cooperation
    S = 0  # Sucker's payoff
    T = 7  # Temptation to defect
    P = 1  # Punishment for mutual defection

    # --- Quantum Game Setup ---
    # For an optimal outcome, we use maximum entanglement.
    gamma = np.pi / 2

    # The initial state for the game is |CC>
    initial_state_vector = np.array([1, 0, 0, 0])

    # The entangling operator J and its conjugate transpose J_dag
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])
    J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(sigma_x, sigma_x)
    J_dag = J.conj().T

    # The available strategies for the players
    # C (Cooperate), D (Defect), and Q (Quantum)
    U_Q = np.array([[1j, 0], [0, -1j]])  # The Quantum strategy that forms the new equilibrium

    # The equilibrium point is where both players choose the Quantum strategy 'Q'
    U_Alice = U_Q
    U_Bob = U_Q

    # --- Calculation Steps ---
    # 1. Entangle the initial state
    psi_entangled = J @ initial_state_vector

    # 2. Players apply their strategies
    U_players = np.kron(U_Alice, U_Bob)
    psi_after_moves = U_players @ psi_entangled

    # 3. Disentangle the state before measurement
    psi_final = J_dag @ psi_after_moves

    # 4. Calculate probabilities of the four classical outcomes
    # The probabilities are the squared magnitudes of the final state vector's components.
    probabilities = np.abs(psi_final)**2
    p_cc, p_cd, p_dc, p_dd = probabilities[0], probabilities[1], probabilities[2], probabilities[3]

    # 5. Calculate the resulting payoff for a player
    player_payoff = R * p_cc + S * p_cd + T * p_dc + P * p_dd

    print("The equilibrium point in the quantum prisoner's dilemma is the (Quantum, Quantum) strategy profile.")
    print("This results in the Pareto optimal outcome, resolving the dilemma.")
    print("The payoff at this equilibrium point is calculated as follows:")
    # The final print statement adheres to the required format of showing the full equation.
    print(f"{R} * {p_cc:.1f} + {S} * {p_cd:.1f} + {T} * {p_dc:.1f} + {P} * {p_dd:.1f} = {player_payoff:.1f}")

solve_quantum_prisoners_dilemma()
<<<5>>>