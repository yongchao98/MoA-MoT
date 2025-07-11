import math

def calculate_growth_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W) 
    growth rates based on true and mistaken probabilities.
    """
    
    # Step 1: Define probabilities and odds
    # True probabilities
    p_true = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Mistaken probabilities
    p_believed = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Net odds (b), where a payout of x:1 means net odds are b = x - 1
    odds = {'A': 3, 'B': 2, 'C': 2}

    # --- Part 1: Calculate the Optimal Growth Rate (W*) ---
    
    # Determine the optimal betting fractions (f_star) using the true probabilities.
    # The Kelly fraction is f = (b*p - (1-p)) / b. We only bet if the edge is positive.
    pA_true, pB_true, pC_true = p_true['A'], p_true['B'], p_true['C']
    bA, bB, bC = odds['A'], odds['B'], odds['C']

    # Check the edge for each bet under true probabilities
    edge_A_true = bA * pA_true - (1 - pA_true) # 3 * 0.5 - 0.5 = 1.0 > 0. Bet.
    edge_B_true = bB * pB_true - (1 - pB_true) # 2 * 0.25 - 0.75 = -0.25. No bet.
    edge_C_true = bC * pC_true - (1 - pC_true) # 2 * 0.25 - 0.75 = -0.25. No bet.

    # Optimal strategy is to only bet on A
    f_star_A = edge_A_true / bA
    f_star_B = 0
    f_star_C = 0

    # Calculate the optimal growth rate W*
    # W* = p_A*log(1 + f_A*b_A) + (p_B + p_C)*log(1 - f_A)
    w_star = pA_true * math.log(1 + f_star_A * bA) + (1 - pA_true) * math.log(1 - f_star_A)

    # --- Part 2: Calculate the Actual Growth Rate (W) ---

    # Determine the bettor's fractions (f) using the mistaken probabilities.
    pA_bel, pB_bel, pC_bel = p_believed['A'], p_believed['B'], p_believed['C']
    
    # Check the edge for each bet under mistaken probabilities
    edge_A_bel = bA * pA_bel - (1 - pA_bel) # 3 * 0.25 - 0.75 = 0. No bet.
    edge_B_bel = bB * pB_bel - (1 - pB_bel) # 2 * 0.5 - 0.5 = 0.5 > 0. Bet.
    edge_C_bel = bC * pC_bel - (1 - pC_bel) # 2 * 0.25 - 0.75 = -0.25. No bet.
    
    # Mistaken strategy is to only bet on B
    f_A = 0
    f_B = edge_B_bel / bB
    f_C = 0

    # Calculate the actual growth rate W using the mistaken bet fractions but true probabilities
    # W = pA_true*log(outcome if A wins) + pB_true*log(outcome if B wins) + pC_true*log(outcome if C wins)
    w = (pA_true * math.log(1 - f_B) +
         pB_true * math.log(1 + f_B * bB) +
         pC_true * math.log(1 - f_B))

    # --- Part 3: Compute the Final Difference ---

    difference = w_star - w

    # The exact symbolic form of the difference is (11/4)*log(2) - (3/2)*log(3)
    print("The final difference can be written with the equation:")
    print("W* - W = (A/B) * log(C) - (D/E) * log(F)\n")
    print("The numbers in this equation are:")
    print("A = 11")
    print("B = 4")
    print("C = 2")
    print("D = 3")
    print("E = 2")
    print("F = 3\n")
    
    print("--- Numerical Results ---")
    print(f"Optimal Growth Rate (W*): {w_star:.4f}")
    print(f"Actual Growth Rate (W):   {w:.4f}")
    print(f"Difference (W* - W):      {difference:.4f}")


calculate_growth_difference()
<<<0.2582>>>