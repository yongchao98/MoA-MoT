import math

def calculate_growth_difference():
    """
    Calculates the difference between optimal (W_star) and actual (W) 
    wealth growth rates based on the Kelly Criterion for a horse race scenario.
    """
    # Step 1: Define probabilities and net odds
    p_true = {'A': 0.5, 'B': 0.25, 'C': 0.25}
    p_mistaken = {'A': 0.25, 'B': 0.5, 'C': 0.25}
    # Payout ratios are 4:1, 3:1, 3:1. Net odds b = payout - 1.
    b = {'A': 3, 'B': 2, 'C': 2}
    competitors = ['A', 'B', 'C']

    # Step 2: Calculate the optimal strategy and W_star
    f_optimal = {c: 0 for c in competitors}
    
    # Using true probabilities, find which bets have an edge (p*(b+1) > 1)
    # Only competitor A has a positive edge: 0.5 * (3 + 1) = 2.0 > 1
    # For B: 0.25 * (2 + 1) = 0.75 < 1
    # For C: 0.25 * (2 + 1) = 0.75 < 1
    # So, we only bet on A.
    p_A = p_true['A']
    b_A = b['A']
    f_optimal['A'] = (p_A * (b_A + 1) - 1) / b_A
    
    # Calculate returns for each outcome with the optimal strategy
    total_f_optimal = sum(f_optimal.values())
    returns_optimal = {}
    for winner in competitors:
        # Wealth factor = (1 - total fraction bet) + fraction_on_winner * (net_odds + 1)
        returns_optimal[winner] = (1 - total_f_optimal) + f_optimal[winner] * (b[winner] + 1)

    # Calculate W_star using true probabilities
    W_star = sum(p_true[c] * math.log(returns_optimal[c]) for c in competitors)

    # Step 3: Calculate the mistaken strategy and W
    f_mistaken = {c: 0 for c in competitors}

    # Using mistaken probabilities, find which bets have an edge
    # For A: 0.25 * (3 + 1) = 1.0. No edge.
    # For B: 0.5 * (2 + 1) = 1.5 > 1. Bet on B.
    # For C: 0.25 * (2 + 1) = 0.75 < 1. No edge.
    p_B_mistaken = p_mistaken['B']
    b_B = b['B']
    f_mistaken['B'] = (p_B_mistaken * (b_B + 1) - 1) / b_B

    # Calculate returns for each outcome with the mistaken strategy
    total_f_mistaken = sum(f_mistaken.values())
    returns_mistaken = {}
    for winner in competitors:
        returns_mistaken[winner] = (1 - total_f_mistaken) + f_mistaken[winner] * (b[winner] + 1)
    
    # Calculate W using true probabilities and mistaken fractions
    W = sum(p_true[c] * math.log(returns_mistaken[c]) for c in competitors)
    
    # Step 4: Output the results
    print("This script calculates the difference between optimal and actual growth rates (W* - W).\n")
    
    # Print calculation for W*
    # W* = pA*ln(R_A) + pB*ln(R_B) + pC*ln(R_C)
    # R_A = 1 + 3*(1/3) - 0 - 0 = 2
    # R_B = 1 - 1/3 = 2/3
    # R_C = 1 - 1/3 = 2/3
    # W* = (1/2)*ln(2) + (1/4)*ln(2/3) + (1/4)*ln(2/3) = (1/2)*ln(4/3)
    print(f"Optimal fraction to bet on A: f_A* = {p_true['A']} - (1-{p_true['A']})/{b['A']} = {f_optimal['A']:.4f}")
    print(f"Optimal growth rate W* = {p_true['A']}*ln({returns_optimal['A']:.4f}) + {p_true['B']}*ln({returns_optimal['B']:.4f}) + {p_true['C']}*ln({returns_optimal['C']:.4f})")
    print(f"W* = {W_star:.5f}\n")
    
    # Print calculation for W
    # W = pA*ln(R'_A) + pB*ln(R'_B) + pC*ln(R'_C)
    # R'_A = 1 - 1/4 = 3/4
    # R'_B = 1 - 0 + 2*(1/4) - 0 = 3/2
    # R'_C = 1 - 1/4 = 3/4
    # W = (1/2)*ln(3/4) + (1/4)*ln(3/2) + (1/4)*ln(3/4)
    print(f"Mistaken fraction to bet on B: f_B' = {p_mistaken['B']} - (1-{p_mistaken['B']})/{b['B']} = {f_mistaken['B']:.4f}")
    print(f"Actual growth rate W = {p_true['A']}*ln({returns_mistaken['A']:.4f}) + {p_true['B']}*ln({returns_mistaken['B']:.4f}) + {p_true['C']}*ln({returns_mistaken['C']:.4f})")
    print(f"W = {W:.5f}\n")
    
    # Print the final difference
    difference = W_star - W
    print(f"The difference is W* - W = {W_star:.5f} - ({W:.5f}) = {difference:.5f}")

calculate_growth_difference()
print(f"\nAnalytically, the result is (11/4)*ln(2) - (3/2)*ln(3).")
print(f"Value: {(11/4)*math.log(2) - (3/2)*math.log(3):.5f}")

# Final Answer
final_answer = (11/4) * math.log(2) - (3/2) * math.log(3)
print(f'<<<{final_answer}>>>')