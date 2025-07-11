from fractions import Fraction

def solve_betting_problem():
    """
    Calculates the achieved growth rate W and the decrease in growth rate Delta W.
    """
    # Step 1: Define probabilities and odds
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q_wrong = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Net odds o_i are (payout - 1)
    # 4-for-1 means payout is 4, net odds are 3.
    net_odds = [3, 2, 6, 6]

    # Step 2: Calculate the optimal strategy and growth rate (W*)
    # Check for positive expectation with true probabilities: p_i * (o_i + 1) > 1
    # Bike 1: (1/2) * (3 + 1) = 2 > 1. Favorable.
    # Bike 2: (1/4) * (2 + 1) = 3/4 < 1. Not favorable.
    # Bike 3: (1/8) * (6 + 1) = 7/8 < 1. Not favorable.
    # Bike 4: (1/8) * (6 + 1) = 7/8 < 1. Not favorable.
    # So, the optimal strategy is to only bet on Bike 1.

    p1 = p_true[0]
    o1 = net_odds[0]
    # Kelly fraction for Bike 1: f1* = p1 - (1-p1)/o1
    f1_star = p1 - (1 - p1) / o1
    # f1_star = 1/2 - (1/2)/3 = 1/2 - 1/6 = 1/3

    # W* = p1*log(1 + f1*o1) + (1-p1)*log(1-f1)
    # W* = (1/2)*log(1 + 1/3*3) + (1/2)*log(1-1/3)
    # W* = (1/2)*log(2) + (1/2)*log(2/3)
    # W* = (1/2)*log(2) + (1/2)*(log(2) - log(3))
    # W* = log(2) - (1/2)*log(3)
    w_star_log2_coeff = Fraction(1)
    w_star_log3_coeff = Fraction(-1, 2)

    # Step 3: Calculate the incorrect strategy and resulting growth rate (W)
    # Check for positive expectation with incorrect probabilities: q_i * (o_i + 1) > 1
    # Bike 1: (1/4) * (3 + 1) = 1. Not favorable.
    # Bike 2: (1/2) * (2 + 1) = 3/2 > 1. Favorable.
    # Bike 3: (1/8) * (6 + 1) = 7/8 < 1. Not favorable.
    # Bike 4: (1/8) * (6 + 1) = 7/8 < 1. Not favorable.
    # So, the incorrect strategy is to only bet on Bike 2.

    q2 = q_wrong[1]
    o2 = net_odds[1]
    # Kelly fraction for Bike 2 based on wrong belief: f2 = q2 - (1-q2)/o2
    f2 = q2 - (1 - q2) / o2
    # f2 = 1/2 - (1/2)/2 = 1/2 - 1/4 = 1/4

    # Calculate achieved growth rate W using f2 but with true probabilities p_i
    # W = p1*log(1-f2) + p2*log(1+f2*o2) + p3*log(1-f2) + p4*log(1-f2)
    # W = (p1+p3+p4)*log(1-f2) + p2*log(1+f2*o2)
    p2 = p_true[1]
    sum_p_others = p_true[0] + p_true[2] + p_true[3] # 1/2 + 1/8 + 1/8 = 3/4
    # W = (3/4)*log(1-1/4) + (1/4)*log(1+1/4*2)
    # W = (3/4)*log(3/4) + (1/4)*log(3/2)
    # W = (3/4)*(log(3)-log(4)) + (1/4)*(log(3)-log(2))
    # W = (3/4)*(log(3)-2*log(2)) + (1/4)*(log(3)-log(2))
    # W = (3/4)*log(3) - (6/4)*log(2) + (1/4)*log(3) - (1/4)*log(2)
    # W = log(3) - (7/4)*log(2)
    w_log2_coeff = Fraction(-7, 4)
    w_log3_coeff = Fraction(1)

    # Step 4: Calculate the decrease in growth rate (Delta W)
    # Delta W = W* - W
    # Delta W = (log(2) - 1/2*log(3)) - (log(3) - 7/4*log(2))
    # Delta W = (1 + 7/4)*log(2) - (1/2 + 1)*log(3)
    # Delta W = (11/4)*log(2) - (3/2)*log(3)
    delta_w_log2_coeff = w_star_log2_coeff - w_log2_coeff
    delta_w_log3_coeff = w_star_log3_coeff - w_log3_coeff

    print("The achieved doubling rate W is given by the equation:")
    print(f"W = ({w_log3_coeff})*log(3) + ({w_log2_coeff})*log(2)")
    print("\nThe decrease in the doubling rate, Delta W = W* - W, is given by the equation:")
    print(f"Delta W = ({delta_w_log2_coeff})*log(2) + ({delta_w_log3_coeff})*log(3)")
    
    # Final answer format
    final_answer = f"{delta_w_log2_coeff}*log(2) - {abs(delta_w_log3_coeff)}*log(3)"
    print(f"\n<<<{final_answer}>>>")

solve_betting_problem()