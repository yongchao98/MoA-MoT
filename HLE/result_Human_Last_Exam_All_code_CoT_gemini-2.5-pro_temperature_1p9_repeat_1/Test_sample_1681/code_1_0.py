import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Solves for the equilibrium point of the Prisoner's Dilemma in a quantum setting.
    """
    # 1. Define the game's parameters
    # Payoff matrix: (Player A's Payoff, Player B's Payoff)
    # (Cooperate, Cooperate), (Cooperate, Defect)
    # (Defect, Cooperate),   (Defect, Defect)
    payoffs = {
        "CC": (5, 5), "CD": (0, 7),
        "DC": (7, 0), "DD": (1, 1)
    }

    # Pauli matrices and Identity
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # 2. Set up the quantum framework (EWL Protocol)
    # The entangling operator J creates a maximally entangled state from |00>
    # J = 1/sqrt(2) * (I ⊗ I + i * σ_x ⊗ σ_x)
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
    J_dag = J.conj().T

    # The initial classical state is |00>, represented as a vector
    psi_00_vector = np.array([1, 0, 0, 0], dtype=complex)

    # The game starts with the entangled state
    psi_entangled = J @ psi_00_vector

    # 3. Define the players' strategies
    # The optimal quantum strategy Q (the "miracle move")
    Q = 1j * sigma_z
    
    # Both players choose the optimal strategy Q
    U_players = np.kron(Q, Q)

    # 4. Calculate the final state and outcome probabilities
    # The final state is calculated by applying the players' strategies
    # and then the disentangling operator J_dag
    psi_final = J_dag @ U_players @ psi_entangled

    # The probabilities of the four classical outcomes are the squared magnitudes
    # of the final state vector's components.
    probabilities = np.abs(psi_final)**2
    prob_cc = probabilities[0]  # Corresponds to |00> (C, C)
    prob_cd = probabilities[1]  # Corresponds to |01> (C, D)
    prob_dc = probabilities[2]  # Corresponds to |10> (D, C)
    prob_dd = probabilities[3]  # Corresponds to |11> (D, D)

    # 5. Calculate the payoffs at equilibrium
    payoff_A = (payoffs["CC"][0] * prob_cc +
                payoffs["CD"][0] * prob_cd +
                payoffs["DC"][0] * prob_dc +
                payoffs["DD"][0] * prob_dd)

    payoff_B = (payoffs["CC"][1] * prob_cc +
                payoffs["CD"][1] * prob_cd +
                payoffs["DC"][1] * prob_dc +
                payoffs["DD"][1] * prob_dd)
    
    # 6. Print the results
    print("In the quantum version of the Prisoner's Dilemma, a new Nash Equilibrium emerges, resolving the classical dilemma.")
    print("This is achieved when both players adopt the quantum strategy 'Q'.")
    print("\n--- Equilibrium Calculation ---")
    print(f"Probabilities of classical outcomes for the (Q, Q) strategy:")
    print(f"  P(Cooperate, Cooperate) = {prob_cc:.4f}")
    print(f"  P(Cooperate, Defect)   = {prob_cd:.4f}")
    print(f"  P(Defect, Cooperate)   = {prob_dc:.4f}")
    print(f"  P(Defect, Defect)      = {prob_dd:.4f}")

    print("\nThis means the (Cooperate, Cooperate) outcome is guaranteed.")
    
    print("\nPlayer A's equilibrium payoff calculation:")
    print(f"Payoff_A = ({payoffs['CC'][0]} * {prob_cc:.2f}) + ({payoffs['CD'][0]} * {prob_cd:.2f}) + ({payoffs['DC'][0]} * {prob_dc:.2f}) + ({payoffs['DD'][0]} * {prob_dd:.2f})")
    print(f"Payoff_A = {payoff_A:.2f}")

    print("\nPlayer B's equilibrium payoff calculation:")
    print(f"Payoff_B = ({payoffs['CC'][1]} * {prob_cc:.2f}) + ({payoffs['CD'][1]} * {prob_cd:.2f}) + ({payoffs['DC'][1]} * {prob_dc:.2f}) + ({payoffs['DD'][1]} * {prob_dd:.2f})")
    print(f"Payoff_B = {payoff_B:.2f}")
    
    print("\nThe equilibrium point is a payoff of (5, 5), which is the Pareto optimal outcome.")

solve_quantum_prisoners_dilemma()