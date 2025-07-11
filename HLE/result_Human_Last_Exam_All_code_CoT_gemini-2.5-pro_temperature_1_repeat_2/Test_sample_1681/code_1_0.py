import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    under maximal entanglement.
    """
    # 1. Define the payoff matrix values from the problem description
    # (C,C) -> (R,R), (C,D) -> (S,T), (D,C) -> (T,S), (D,D) -> (P,P)
    R = 5  # Reward for mutual cooperation
    S = 0  # Sucker's payoff for cooperating while opponent defects
    T = 7  # Temptation payoff for defecting while opponent cooperates
    P = 1  # Punishment for mutual defection

    # 2. Define the quantum states and operators
    # Single qubit basis states |C> = |0>, |D> = |1>
    C = np.array([[1], [0]])
    D = np.array([[0], [1]])

    # Initial two-qubit state is |CC>
    initial_state_CC = np.kron(C, C)

    # Standard operators
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])

    # 3. Set up the quantum game protocol (EWL)
    # We choose maximal entanglement for an optimal outcome, as requested.
    gamma = np.pi / 2

    # Entangling operator J(gamma)
    J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(sigma_x, sigma_x)
    # Its conjugate transpose (the disentangling operator)
    J_dag = J.conj().T

    # 4. Define the equilibrium strategies
    # The new Nash Equilibrium strategy in the quantum game is 'Q'
    Q_strategy = np.array([[1j, 0], [0, -1j]])

    # Both players adopt the optimal quantum strategy
    U_A = Q_strategy
    U_B = Q_strategy
    
    # The full operator for the players' moves
    players_operator = np.kron(U_A, U_B)

    # 5. Simulate the game's evolution
    # a. The initial state is entangled: |ψ_i> = J|CC>
    psi_initial = J @ initial_state_CC
    
    # b. Players apply their strategies: |ψ_moves> = (U_A ⊗ U_B)|ψ_i>
    psi_after_moves = players_operator @ psi_initial

    # c. The state is disentangled before measurement: |ψ_final> = J†|ψ_moves>
    psi_final = J_dag @ psi_after_moves

    # 6. Calculate outcome probabilities from the final state vector
    # Probabilities are the squared magnitudes of the final state's amplitudes
    # psi_final = [α_cc, α_cd, α_dc, α_dd]^T
    prob_CC = np.abs(psi_final[0, 0])**2
    prob_CD = np.abs(psi_final[1, 0])**2
    prob_DC = np.abs(psi_final[2, 0])**2
    prob_DD = np.abs(psi_final[3, 0])**2

    # 7. Calculate the final payoffs
    # Since the quantum result is exact, we can round to get clean integers for the equation.
    p_cc, p_cd, p_dc, p_dd = [int(round(p)) for p in [prob_CC, prob_CD, prob_DC, prob_DD]]

    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd

    # 8. Print the results including the full payoff equations
    print("Under optimal entanglement, a new Nash Equilibrium (Q, Q) emerges.")
    print("The final probabilities for the classical outcomes are:")
    print(f"P(Cooperate, Cooperate) = {p_cc}")
    print(f"P(Cooperate, Defect)  = {p_cd}")
    print(f"P(Defect,   Cooperate) = {p_dc}")
    print(f"P(Defect,   Defect)  = {p_dd}\n")

    print("The payoff for each player is calculated as follows:")
    print(f"Player A's Payoff: {R} * {p_cc} + {S} * {p_cd} + {T} * {p_dc} + {P} * {p_dd} = {payoff_A}")
    print(f"Player B's Payoff: {R} * {p_cc} + {T} * {p_cd} + {S} * {p_dc} + {P} * {p_dd} = {payoff_B}\n")
    
    print(f"The equilibrium point is the payoff pair: ({payoff_A}, {payoff_B})")


if __name__ == "__main__":
    solve_quantum_prisoners_dilemma()

<<<(5, 5)>>>