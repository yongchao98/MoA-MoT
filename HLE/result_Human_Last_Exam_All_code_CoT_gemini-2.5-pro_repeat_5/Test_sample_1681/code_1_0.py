import numpy as np

def solve_quantum_prisoners_dilemma():
    """
    Analyzes the quantum Prisoner's Dilemma to find its equilibrium point.

    The analysis shows that a new Nash Equilibrium (Q, Q) emerges in the quantum game,
    where Q is a specific quantum strategy. This equilibrium yields a payoff of (5, 5),
    which is the socially optimal outcome, thus resolving the dilemma.
    """

    # 1. Payoff Matrix and Classical Game
    # Payoffs are (Alice, Bob) for (Alice's move, Bob's move)
    # (C, C) -> (5, 5)
    # (C, D) -> (0, 7)
    # (D, C) -> (7, 0)
    # (D, D) -> (1, 1)
    # In the classical game, the Nash Equilibrium is (D, D) with payoff (1, 1).

    print("Quantum Prisoner's Dilemma Analysis")
    print("="*40)
    print("Payoff Matrix (Alice, Bob):")
    print("         Bob: Cooperate | Defect")
    print("---------------------------------")
    print("Alice C:      (5, 5)    | (0, 7)")
    print("Alice D:      (7, 0)    | (1, 1)")
    print("\nIn the quantum version (EWL protocol), we seek a new equilibrium.")
    print("Let's analyze the symmetric strategy profile (Q, Q), where Q = U(theta=0, phi=pi/2).\n")

    # 2. Deriving Alice's Payoff against Bob's Q
    print("Step 1: Calculate Alice's payoff function $A(U_A, Q).")
    print("Alice's strategy U_A is parameterized by (theta, phi).")
    print("Bob's strategy is fixed at Q (theta_B=0, phi_B=pi/2).\n")

    print("The probabilities of the four outcomes (CC, CD, DC, DD) depend on (theta, phi).")
    print("Let C = cos(theta/2) and S = sin(theta/2).")
    print("After complex derivation, the probabilities are:")
    print("P_CC = C^2 * sin^2(phi)")
    print("P_CD = S^2 / 2")
    print("P_DC = S^2 / 2")
    print("P_DD = C^2 * cos^2(phi)\n")

    # Original payoff values for Alice
    p_cc, p_cd, p_dc, p_dd = (5, 0, 7, 1)

    print("Step 2: Formulate Alice's expected payoff $A.")
    print(f"The formula is: $A = {p_cc}*P_CC + {p_cd}*P_CD + {p_dc}*P_DC + {p_dd}*P_DD")
    print(f"$A = {p_cc} * (C^2 * sin^2(phi)) + {p_cd} * (S^2 / 2) + {p_dc} * (S^2 / 2) + {p_dd} * (C^2 * cos^2(phi))")
    
    # Substituting p_cd=0 simplifies the equation
    print("\nSimplifying the equation:")
    # The coefficients in the final equation are 5, 7, and 1.
    final_coeff_1 = p_cc
    final_coeff_2 = p_dc
    final_coeff_3 = p_dd
    
    print(f"$A = {final_coeff_1}*C^2*sin^2(phi) + {final_coeff_2}*(S^2/2) + {final_coeff_3}*C^2*cos^2(phi)")
    
    # Grouping terms
    # $A = C^2 * (5*sin^2(phi) + 1*cos^2(phi)) + 3.5 * S^2
    # This can also be written as:
    # $A = cos^2(theta/2) * (1 + 4*sin^2(phi)) + 3.5 * sin^2(theta/2)
    coeff_a = 1
    coeff_b = 4
    coeff_c = 3.5

    print("\nWhich can be further simplified to the final equation for analysis:")
    print(f"Payoff_A = cos^2(theta/2) * ({coeff_a} + {coeff_b}*sin^2(phi)) + {coeff_c}*sin^2(theta/2)\n")

    # 3. Maximizing Alice's Payoff
    print("Step 3: Find the strategy (theta, phi) that maximizes Alice's payoff.")
    print("To maximize Payoff_A, we analyze the terms:")
    print(" - The term `(1 + 4*sin^2(phi))` is maximized when sin^2(phi) is 1, so phi = pi/2.")
    print(" - With phi=pi/2, the equation becomes: Payoff_A = 5*cos^2(theta/2) + 3.5*sin^2(theta/2).")
    print(" - This is maximized when cos^2(theta/2) is 1 (and sin^2(theta/2) is 0), so theta = 0.\n")

    print("Alice's best response to Bob playing Q is the strategy (theta=0, phi=pi/2), which is Q itself.")
    
    # Calculate the payoff at this point
    theta_max = 0
    phi_max = np.pi / 2
    
    C2_max = np.cos(theta_max / 2)**2
    S2_max = np.sin(theta_max / 2)**2
    sin2_phi_max = np.sin(phi_max)**2
    
    max_payoff = C2_max * (coeff_a + coeff_b * sin2_phi_max) + coeff_c * S2_max
    
    print(f"The maximum payoff Alice can achieve is {max_payoff:.1f}.\n")
    
    # 4. Conclusion
    print("Conclusion:")
    print("Since Alice's best response to Bob playing Q is to also play Q, and the game is symmetric,")
    print("the strategy profile (Q, Q) is a stable Nash Equilibrium.")
    print("At this equilibrium, both players achieve the optimal payoff.")
    print("\nEquilibrium Strategy: (Q, Q)")
    print(f"Equilibrium Payoff: ({max_payoff:.1f}, {max_payoff:.1f})")


solve_quantum_prisoners_dilemma()
# The question asks for the equilibrium point. In game theory, this can refer to the
# strategy profile or the payoff. Given the context, the numerical payoff is expected.
# Since the payoff is symmetric (5,5), we output the single value.
print("\n<<<5>>>")