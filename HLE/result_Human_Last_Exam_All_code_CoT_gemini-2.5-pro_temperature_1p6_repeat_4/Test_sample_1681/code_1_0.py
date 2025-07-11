import numpy as np

def print_payoff_calculation(player_name, payoffs, probs):
    """Prints the detailed payoff calculation."""
    p_cc, p_cd, p_dc, p_dd = probs
    r, s, t, p = payoffs
    total_payoff = p_cc * r + p_cd * s + p_dc * t + p_dd * p
    print(f"  {player_name}'s payoff equation: {r:.1f}*{p_cc:.2f} + {s:.1f}*{p_cd:.2f} + {t:.1f}*{p_dc:.2f} + {p:.1f}*{p_dd:.2f} = {total_payoff:.2f}")
    return total_payoff

def calculate_payoffs(Ua, Ub, psi_in, J_dag, payoff_matrix):
    """Calculates the final state and expected payoffs for players A and B."""
    
    # Combine players' strategies
    U_combined = np.kron(Ua, Ub)
    
    # Calculate the final state before measurement
    psi_f = J_dag @ U_combined @ psi_in
    
    # Calculate probabilities of the four outcomes |CC>, |CD>, |DC>, |DD>
    probs = np.abs(psi_f.flatten())**2
    p_cc, p_cd, p_dc, p_dd = probs
    
    # Payoffs from the matrix: R, S, T, P
    # [[(R,R) , (S,T)],
    #  [(T,S) , (P,P)]]
    R = payoff_matrix[0, 0, 0] # Cooperate-Cooperate (Reward)
    S = payoff_matrix[0, 1, 0] # Cooperate-Defect (Sucker)
    T = payoff_matrix[1, 0, 0] # Defect-Cooperate (Temptation)
    P = payoff_matrix[1, 1, 0] # Defect-Defect (Punishment)

    payoffs_A = np.array([R, S, T, P])
    payoffs_B = np.array([R, T, S, P])

    print(f"  Probabilities (CC, CD, DC, DD): ({p_cc:.2f}, {p_cd:.2f}, {p_dc:.2f}, {p_dd:.2f})")
    
    payoff_A = print_payoff_calculation("Player A", payoffs_A, probs)
    payoff_B = print_payoff_calculation("Player B", payoffs_B, probs)
    
    return (payoff_A, payoff_B)


def solve_quantum_dilemma():
    """
    Sets up and solves the Quantum Prisoner's Dilemma, printing the analysis.
    """
    # Payoff Matrix [[(A,B), (A,B)], [(A,B), (A,B)]]
    # Upper-left: (Cooperate, Cooperate), Bottom-right: (Defect, Defect)
    payoff_matrix = np.array([
        [(5, 5), (0, 7)],
        [(7, 0), (1, 1)]
    ])

    # Define quantum operators (Pauli matrices and Identity)
    I = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)

    # Define player strategies
    C = I              # Cooperate
    D = 1j * sx        # Defect
    Q = 1j * sz        # Quantum
    
    strategies = {'C': C, 'D': D, 'Q': Q}

    # Define initial state and entangling/disentangling operators
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sx, sx))
    J_dag = J.conj().T
    psi_00 = np.array([1, 0, 0, 0], dtype=complex).reshape(4, 1)
    psi_in = J @ psi_00

    print("Analyzing the Quantum Prisoner's Dilemma:\n")

    # 1. Classical Dilemma Outcome
    print("1. Classical Equilibrium (D, D):")
    p_dd = calculate_payoffs(strategies['D'], strategies['D'], psi_in, J_dag, payoff_matrix)
    print(f"-> The classical equilibrium (Defect, Defect) yields a payoff of ({p_dd[0]:.1f}, {p_dd[1]:.1f}).\n")

    # 2. Proposed Quantum Solution
    print("2. Quantum Strategy Equilibrium (Q, Q):")
    p_qq = calculate_payoffs(strategies['Q'], strategies['Q'], psi_in, J_dag, payoff_matrix)
    print(f"-> The quantum strategy (Q, Q) yields a superior payoff of ({p_qq[0]:.1f}, {p_qq[1]:.1f}).\n")

    # 3. Test if (Q, Q) is a Nash Equilibrium
    print("3. Checking for profitable deviations from (Q, Q):")
    print("   If Player B plays Q, does Player A have an incentive to switch?")
    
    # A switches to C
    print("\n   Case A: Player A switches from Q to C (Game is C, Q)")
    p_cq = calculate_payoffs(strategies['C'], strategies['Q'], psi_in, J_dag, payoff_matrix)
    print(f"   -> Player A's payoff would change from {p_qq[0]:.1f} to {p_cq[0]:.1f}. This is not a better outcome for A.")

    # A switches to D
    print("\n   Case B: Player A switches from Q to D (Game is D, Q)")
    p_dq = calculate_payoffs(strategies['D'], strategies['Q'], psi_in, J_dag, payoff_matrix)
    print(f"   -> Player A's payoff would change from {p_qq[0]:.1f} to {p_dq[0]:.1f}. This is not a better outcome for A.")

    print("\nConclusion:")
    print("Neither player has an incentive to unilaterally deviate from strategy Q if the other is playing Q.")
    print("Therefore, (Q, Q) is a Nash Equilibrium.")
    print("The equilibrium point corresponds to the strategy pair (Q, Q), which yields the optimal payoff of (5, 5).")

if __name__ == '__main__':
    solve_quantum_dilemma()