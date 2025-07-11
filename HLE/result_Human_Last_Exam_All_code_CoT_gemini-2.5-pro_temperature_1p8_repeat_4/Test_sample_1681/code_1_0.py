import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Solves for the equilibrium point of a quantum Prisoner's Dilemma.

    The classical payoff matrix is:
    [[(5,5) , (0,7)],  (C,C), (C,D)
     [(7,0) , (1,1)]]  (D,C), (D,D)

    When played in a quantum space with maximal entanglement (using the EWL protocol),
    the game transforms into a different game. For the strategies Cooperate (C) and Defect (D),
    the new payoff matrix has the structure of a "Game of Chicken".
    """
    print("Step 1: Define the payoffs for the transformed quantum game (Game of Chicken).")
    # These payoffs are derived from the quantum game theory formalism with gamma = pi/2.
    # P_alicestrat_bobstrat[player_index]
    P_CC = (3.0, 3.0)
    P_CD = (3.5, 3.5)
    P_DC = (3.5, 3.5)
    P_DD = (3.0, 3.0)

    print(f"Payoff(C,C) = {P_CC}")
    print(f"Payoff(C,D) = {P_CD}")
    print(f"Payoff(D,C) = {P_DC}")
    print(f"Payoff(D,D) = {P_DD}\n")

    print("Step 2: Find the mixed strategy Nash Equilibrium.")
    print("Let p be the probability Alice plays Defect, and (1-p) she plays Cooperate.")
    print("Let q be the probability Bob plays Defect, and (1-q) he plays Cooperate.")
    print("For a mixed NE, each player must be indifferent between their pure strategies.\n")
    
    print("Alice's perspective: Expected Payoff(playing C) = Expected Payoff(playing D)")
    # Payoff_Alice(C) = q * P_CD[0] + (1-q) * P_CC[0]
    # Payoff_Alice(D) = q * P_DD[0] + (1-q) * P_DC[0]
    print(f"q * {P_CD[0]} + (1-q) * {P_CC[0]} = q * {P_DD[0]} + (1-q) * {P_DC[0]}")
    # We solve for q
    # P_CC[0] - P_DC[0] = q * (P_DD[0] - P_DC[0] - P_CD[0] + P_CC[0])
    q_numerator = P_CC[0] - P_DC[0]
    q_denominator = P_DD[0] - P_DC[0] - P_CD[0] + P_CC[0]
    q = q_numerator / q_denominator
    print(f"Solving for q: q = ({P_CC[0]} - {P_DC[0]}) / ({P_DD[0]} - {P_DC[0]} - {P_CD[0]} + {P_CC[0]}) = {q_numerator} / {q_denominator} = {q}")
    # Due to symmetry, p = q
    p = q
    print(f"By symmetry, p = {p}\n")
    
    print("Step 3: Calculate the expected payoff at this equilibrium.")
    # Expected Payoff = Sum of (Prob(outcome) * Payoff(outcome))
    # E = p*q*P_DD + p*(1-q)*P_DC + (1-p)*q*P_CD + (1-p)*(1-q)*P_CC
    p_d_d = p * q
    p_d_c = p * (1-q)
    p_c_d = (1-p) * q
    p_c_c = (1-p) * (1-q)

    expected_payoff = (p_c_c * P_CC[0] +
                       p_c_d * P_CD[0] +
                       p_d_c * P_DC[0] +
                       p_d_d * P_DD[0])

    print("The final equilibrium payoff is the expected value of the game for a player.")
    print("E[$] = P(C,C)*Payoff(C,C) + P(C,D)*Payoff(C,D) + P(D,C)*Payoff(D,C) + P(D,D)*Payoff(D,D)")
    print("P(C,C) = (1-p)*(1-q) =", p_c_c)
    print("P(C,D) = (1-p)*q =", p_c_d)
    print("P(D,C) = p*(1-q) =", p_d_c)
    print("P(D,D) = p*q =", p_d_d)
    print(f"E[$] = ({p_c_c:.2f} * {P_CC[0]}) + ({p_c_d:.2f} * {P_CD[0]}) + ({p_d_c:.2f} * {P_DC[0]}) + ({p_d_d:.2f} * {P_DD[0]}) = {expected_payoff:.4f}")

    print("\nThe equilibrium point (expected payoff) is:")
    print(expected_payoff)
    return expected_payoff

final_payoff = solve_quantum_prisoners_dilemma()
# The required output format.
# <<<answer content>>>
# For example: <<<C>>>, <<<9.8>>>
print(f'<<<{final_payoff}>>>')
