import numpy as np

def solve_quantum_dilemma():
    """
    Analyzes the Quantum Prisoner's Dilemma to find the equilibrium point.
    """
    # Step 1: Define the payoff matrix from the problem description
    #         Player 2
    #         Cooperate (C)   Defect (D)
    # P1 C      (R, R)          (S, T)
    #    D      (T, S)          (P, P)
    R = 5  # Reward for mutual cooperation
    S = 0  # Sucker's payoff
    T = 7  # Temptation to defect
    P = 1  # Punishment for mutual defection
    payoffs = {'R': R, 'S': S, 'T': T, 'P': P}

    # Step 2: Define the quantum mechanical setup
    # Basis vectors for the two-qubit system
    ket00 = np.array([[1], [0], [0], [0]], dtype=complex)

    # Entangling operator J for maximal entanglement
    I2 = np.identity(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    J = (1 / np.sqrt(2)) * (np.kron(I2, I2) + 1j * np.kron(sx, sx))

    # Step 3: Define player strategies as unitary operators
    # The 'miracle' move that allows for cooperation
    Q = np.array([[1j, 0], [0, -1j]], dtype=complex)

    # The optimal counter-strategy to Q, a quantum form of 'Defect'
    D_prime = np.array([[0, 1], [-1, 0]], dtype=complex)

    # Step 4: Analyze the game to find the Nash Equilibrium
    # The game played with strategies Q and D' results in a classic Prisoner's Dilemma,
    # where (D', D') is the dominant strategy and the only Nash Equilibrium.
    # We will calculate the payoff for this equilibrium point.
    equilibrium_strategy_A = D_prime
    equilibrium_strategy_B = D_prime

    # Step 5: Calculate the final state and payoffs for the equilibrium
    # The final state |ψ'> is calculated as J† (U_A ⊗ U_B) J |ψ_initial>
    psi_initial = J @ ket00
    U_game = np.kron(equilibrium_strategy_A, equilibrium_strategy_B)
    psi_transformed = U_game @ psi_initial
    J_dag = J.conj().T
    psi_final = J_dag @ psi_transformed

    # Probabilities of the four classical outcomes (CC, CD, DC, DD)
    prob_CC = np.abs(psi_final[0, 0])**2
    prob_CD = np.abs(psi_final[1, 0])**2
    prob_DC = np.abs(psi_final[2, 0])**2
    prob_DD = np.abs(psi_final[3, 0])**2
    
    # Round for clean display
    (pCC, pCD, pDC, pDD) = (round(p.real, 2) for p in [prob_CC, prob_CD, prob_DC, prob_DD])

    # Calculate expected payoffs for each player
    payoff_A = R * pCC + S * pCD + T * pDC + P * pDD
    payoff_B = R * pCC + T * pCD + S * pDC + P * pDD

    # Step 6: Print the detailed calculation for the equilibrium payoff
    print("In the quantum Prisoner's Dilemma, the Nash Equilibrium strategy is for both players to play a defecting strategy D'.")
    print("Here is the calculation for the payoff at the equilibrium point (D', D'):\n")
    print(f"Probabilities of outcomes: P(Cooperate, Cooperate) = {pCC}, P(Cooperate, Defect) = {pCD}, P(Defect, Cooperate) = {pDC}, P(Defect, Defect) = {pDD}\n")
    
    print("Final Payoff Equation for Player A:")
    print(f"Payoff(A) = (R * P(CC)) + (S * P(CD)) + (T * P(DC)) + (P * P(DD))")
    print(f"          = ({R} * {pCC}) + ({S} * {pCD}) + ({T} * {pDC}) + ({P} * {pDD})")
    print(f"          = {payoff_A:.2f}\n")

    print("Final Payoff Equation for Player B:")
    print(f"Payoff(B) = (R * P(CC)) + (T * P(CD)) + (S * P(DC)) + (P * P(DD))")
    print(f"          = ({R} * {pCC}) + ({T} * {pCD}) + ({S} * {pDC}) + ({P} * {pDD})")
    print(f"          = {payoff_B:.2f}\n")
    
    print(f"The equilibrium point is the payoff profile ({payoff_A:.0f}, {payoff_B:.0f}).")

solve_quantum_dilemma()