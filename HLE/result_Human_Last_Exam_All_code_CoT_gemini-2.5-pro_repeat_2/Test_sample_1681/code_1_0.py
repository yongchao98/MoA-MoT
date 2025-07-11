import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates and identifies the Nash Equilibrium for the Quantum Prisoner's Dilemma.
    """
    # Payoff matrix from the problem description
    #           B: Cooperate(0)  B: Defect(1)
    # A: Coop(0)   (5, 5)           (0, 7)
    # A: Defect(1)  (7, 0)           (1, 1)
    payoffs_A = np.array([[5, 0], [7, 1]])
    payoffs_B = np.array([[5, 7], [0, 1]])

    # Initial state |00>
    psi_0 = np.array([1, 0, 0, 0], dtype=complex)

    # Maximal entanglement operator J(gamma) for gamma = pi/2
    gamma = np.pi / 2
    J = np.cos(gamma / 2) * np.identity(4) + 1j * np.sin(gamma / 2) * np.kron(
        np.array([[0, 1], [1, 0]]),  # sigma_x
        np.array([[0, 1], [1, 0]])   # sigma_x
    )
    J_dag = J.conj().T

    # Define strategy operators
    # C: Cooperate (Identity)
    U_C = np.identity(2, dtype=complex)
    # D: Defect (a specific unitary operator, U(pi, 0))
    U_D = np.array([[0, 1], [-1, 0]], dtype=complex)
    # Q: Quantum "Miracle" Move (U(0, pi/2))
    U_Q = np.array([[1j, 0], [0, -1j]], dtype=complex)

    strategies = {'C': U_C, 'D': U_D, 'Q': U_Q}
    strategy_names = ['C', 'D', 'Q']
    
    print("--- Quantum Prisoner's Dilemma Payoff Matrix ---")
    print("Payoffs are shown as (Player A, Player B)\n")
    
    # Store results for analysis
    results = {}

    for name_A in strategy_names:
        for name_B in strategy_names:
            U_A = strategies[name_A]
            U_B = strategies[name_B]

            # Combined operator for players' strategies
            U = np.kron(U_A, U_B)

            # Calculate the final state before measurement
            psi_final = J_dag @ U @ J @ psi_0

            # Probabilities of the four classical outcomes |00>, |01>, |10>, |11>
            probs = np.abs(psi_final)**2
            p_cc, p_cd, p_dc, p_dd = probs[0], probs[1], probs[2], probs[3]

            # Calculate expected payoffs
            payoff_A = p_cc * payoffs_A[0, 0] + p_cd * payoffs_A[0, 1] + \
                       p_dc * payoffs_A[1, 0] + p_dd * payoffs_A[1, 1]
            payoff_B = p_cc * payoffs_B[0, 0] + p_cd * payoffs_B[0, 1] + \
                       p_dc * payoffs_B[1, 0] + p_dd * payoffs_B[1, 1]
            
            results[(name_A, name_B)] = (round(payoff_A, 2), round(payoff_B, 2))

    # Print the payoff table
    header = "B->".ljust(6) + "".join([f"{name.center(12)}" for name in strategy_names])
    print(header)
    print("-" * len(header))
    for name_A in strategy_names:
        row_str = f"A:{name_A}".ljust(6)
        for name_B in strategy_names:
            row_str += f"{str(results[(name_A, name_B)]).center(12)}"
        print(row_str)
    
    print("\n--- Analysis ---")
    print("In the classical game, the Nash Equilibrium is (D, D) with payoff (1, 1).")
    print("In this quantum game with maximal entanglement:")
    print("1. The classical (D, D) equilibrium is no longer stable. If player A plays D, player B's best move is Q (payoff 3.5 vs 1).")
    print("2. A new Nash Equilibrium emerges at (Q, Q). If player A plays Q, player B's best move is also Q (payoff 5 vs 3.5 for D or 1 for C).")
    print("This new equilibrium (Q, Q) resolves the dilemma, leading to a Pareto-optimal outcome where both players cooperate via quantum strategies.")
    
    print("\n--- Equilibrium Point Calculation ---")
    print("The equilibrium point is the payoff at the (Q, Q) Nash Equilibrium.")
    
    # Re-calculate for Q,Q to show the equation
    U_A = strategies['Q']
    U_B = strategies['Q']
    U = np.kron(U_A, U_B)
    psi_final = J_dag @ U @ J @ psi_0
    probs = np.abs(psi_final)**2
    p_cc, p_cd, p_dc, p_dd = probs[0], probs[1], probs[2], probs[3]
    
    payoff_A = p_cc * payoffs_A[0, 0] + p_cd * payoffs_A[0, 1] + \
               p_dc * payoffs_A[1, 0] + p_dd * payoffs_A[1, 1]
    payoff_B = p_cc * payoffs_B[0, 0] + p_cd * payoffs_B[0, 1] + \
               p_dc * payoffs_B[1, 0] + p_dd * payoffs_B[1, 1]
               
    print("\nFor the (Q, Q) strategy pair, the outcome probabilities are:")
    print(f"P(Cooperate, Cooperate) = {p_cc:.2f}")
    print(f"P(Cooperate, Defect)   = {p_cd:.2f}")
    print(f"P(Defect,   Cooperate) = {p_dc:.2f}")
    print(f"P(Defect,   Defect)   = {p_dd:.2f}")

    print("\nPlayer A Final Payoff Equation:")
    print(f"{payoff_A:.1f} = {payoffs_A[0,0]}*{p_cc:.2f} + {payoffs_A[0,1]}*{p_cd:.2f} + {payoffs_A[1,0]}*{p_dc:.2f} + {payoffs_A[1,1]}*{p_dd:.2f}")

    print("\nPlayer B Final Payoff Equation:")
    print(f"{payoff_B:.1f} = {payoffs_B[0,0]}*{p_cc:.2f} + {payoffs_B[0,1]}*{p_cd:.2f} + {payoffs_B[1,0]}*{p_dc:.2f} + {payoffs_B[1,1]}*{p_dd:.2f}")

solve_quantum_prisoners_dilemma()