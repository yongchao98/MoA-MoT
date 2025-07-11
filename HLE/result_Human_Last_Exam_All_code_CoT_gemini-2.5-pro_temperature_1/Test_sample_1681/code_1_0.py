import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Calculates the mixed strategy Nash Equilibrium payoff for the Quantum Prisoner's Dilemma
    using the EWL protocol with maximal entanglement.
    """
    # 1. Payoff matrix from the problem description
    # (Cooperate, Cooperate), (Cooperate, Defect)
    # (Defect, Cooperate)  , (Defect, Defect)
    # R=Reward, S=Sucker, T=Temptation, P=Punishment
    R, S, T, P = 5, 0, 7, 1

    # 2. Set up the quantum game parameters (EWL protocol)
    # gamma is the entanglement parameter. We choose maximal entanglement.
    gamma = np.pi / 2
    c = np.cos(gamma / 2)
    s = np.sin(gamma / 2)
    c2 = c**2
    s2 = s**2

    # 3. Calculate the payoffs for the 2x2 quantum meta-game between C and D strategies.
    # These formulas are derived from the EWL protocol.
    # Payoff for Alice when the outcome is (Alice's Strategy, Bob's Strategy)
    
    # Outcome (C, C): P_CC=c^2, P_CD=0, P_DC=0, P_DD=s^2
    payoff_A_CC = R * c2 + S * 0 + T * 0 + P * s2
    
    # Outcome (C, D): P_CC=0, P_CD=c^2, P_DC=s^2, P_DD=0
    payoff_A_CD = R * 0 + S * c2 + T * s2 + P * 0
    
    # Outcome (D, C): P_CC=0, P_CD=s^2, P_DC=c^2, P_DD=0
    payoff_A_DC = R * 0 + S * s2 + T * c2 + P * 0
    
    # Outcome (D, D): P_CC=s^2, P_CD=0, P_DC=0, P_DD=c^2
    payoff_A_DD = R * s2 + S * 0 + T * 0 + P * c2

    # The resulting 2x2 payoff matrix for Alice is:
    # [[payoff_A_CC, payoff_A_CD],
    #  [payoff_A_DC, payoff_A_DD]]

    # 4. Solve for the mixed strategy Nash Equilibrium probability 'p' for playing 'C'.
    # In a symmetric 2x2 game, the mixed NE probability 'p' for player 1 (Alice)
    # is the one that makes player 2 (Bob) indifferent between their two strategies.
    # Due to symmetry, p_Alice = p_Bob.
    # p = (payoff_A_DD - payoff_A_CD) / (payoff_A_CC - payoff_A_DC - payoff_A_CD + payoff_A_DD)
    numerator = payoff_A_DD - payoff_A_CD
    denominator = (payoff_A_CC - payoff_A_DC) - (payoff_A_CD - payoff_A_DD)
    
    if denominator == 0:
        print("No unique mixed strategy equilibrium.")
        return

    p_cooperate = numerator / denominator

    # 5. Calculate the expected payoff for a player at this equilibrium.
    # Expected Payoff = p * (Payoff if opponent plays C) + (1-p) * (Payoff if opponent plays D)
    # We can calculate this from Alice's perspective, assuming Bob also plays C with probability p.
    # E = p * payoff_A_CC + (1-p) * payoff_A_CD
    equilibrium_payoff = p_cooperate * payoff_A_CC + (1 - p_cooperate) * payoff_A_CD

    # The problem asks to output the final equation with numbers.
    print("The game between quantum strategies C and D has the following payoff matrix for Player 1 (Alice):")
    print(f"          Bob: C              Bob: D")
    print(f"Alice: C  ({payoff_A_CC:.2f}, ...)        ({payoff_A_CD:.2f}, ...)")
    print(f"Alice: D  ({payoff_A_DC:.2f}, ...)        ({payoff_A_DD:.2f}, ...)")
    print("\nThis game has a mixed strategy Nash Equilibrium.")
    print(f"The equilibrium probability of choosing 'Cooperate' is p = {p_cooperate:.2f}.")
    print("\nThe equilibrium point (expected payoff) is calculated as follows:")
    print(f"Expected Payoff = p * Payoff(Alice plays C | Bob plays C) + (1-p) * Payoff(Alice plays C | Bob plays D)")
    print(f"Expected Payoff = {p_cooperate:.2f} * {payoff_A_CC:.2f} + (1 - {p_cooperate:.2f}) * {payoff_A_CD:.2f}")
    final_calc_str = f"= {p_cooperate * payoff_A_CC:.3f} + { (1 - p_cooperate) * payoff_A_CD:.3f}"
    print(f"                 {final_calc_str}")
    print(f"                 = {equilibrium_payoff:.4f}")
    
    # Return the final numeric answer as requested by the format.
    return equilibrium_payoff

# Run the calculation and store the result
final_payoff = solve_quantum_prisoners_dilemma()
print(f"\n<<<>>>\n{final_payoff}\n<<<>>>")
