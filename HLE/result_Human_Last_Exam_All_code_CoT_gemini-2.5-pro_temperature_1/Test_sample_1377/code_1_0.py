import math
from fractions import Fraction

def solve_growth_rate_difference():
    """
    Calculates the difference between the optimal (W*) and actual (W)
    wealth growth rates based on the Kelly Criterion.
    """
    # Step 1: Define problem parameters
    # True probabilities for (A, B, C)
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 4)]
    # Mistakenly believed probabilities for (A, B, C)
    p_believed = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)]
    # Net odds b = (payout - 1) for (A, B, C)
    b = [3, 2, 2]

    # Step 2: Calculate the optimal strategy and growth rate (W*)
    # We only bet if the edge p*b - (1-p) is positive.
    # Using true probabilities:
    # Edge A: (1/2)*3 - (1/2) = 1 > 0. Bet.
    # Edge B: (1/4)*2 - (3/4) = -1/4 < 0. No bet.
    # Edge C: (1/4)*2 - (3/4) = -1/4 < 0. No bet.
    # The optimal strategy is to only bet on A.
    
    # Kelly fraction f = (p*b - q) / b, where q = 1-p
    p_A_true = p_true[0]
    b_A = b[0]
    f_star_A = (p_A_true * b_A - (1 - p_A_true)) / b_A

    # W* is the expected log return from the optimal strategy
    # If A wins (prob 1/2), return is 1 + f*b = 1 + (1/3)*3 = 2
    # If A loses (prob 1/2), return is 1 - f = 1 - 1/3 = 2/3
    W_star = p_A_true * math.log(1 + f_star_A * b_A) + (1 - p_A_true) * math.log(1 - f_star_A)

    print("Step 1: Calculate the optimal growth rate W*")
    print(f"The optimal strategy is to bet f_A = {f_star_A:.3f} (or {Fraction(f_star_A).limit_denominator()}) on Competitor A.")
    print(f"The formula for W* is: P(A_wins) * log(1 + f_A*b_A) + P(A_loses) * log(1 - f_A)")
    print(f"W* = {p_A_true} * log(1 + {Fraction(f_star_A).limit_denominator()} * {b_A}) + {1-p_A_true} * log(1 - {Fraction(f_star_A).limit_denominator()})")
    print(f"W* = {p_A_true} * log({1 + f_star_A * b_A}) + {1-p_A_true} * log({1 - f_star_A})")
    print(f"W* = {W_star:.5f}\n")
    
    # Step 3: Calculate the mistaken strategy and actual growth rate (W)
    # Using believed probabilities:
    # Edge A: (1/4)*3 - (3/4) = 0. No bet.
    # Edge B: (1/2)*2 - (1/2) = 1/2 > 0. Bet.
    # Edge C: (1/4)*2 - (3/4) = -1/4 < 0. No bet.
    # The mistaken strategy is to only bet on B.
    
    p_B_believed = p_believed[1]
    b_B = b[1]
    f_mistaken_B = (p_B_believed * b_B - (1 - p_B_believed)) / b_B
    
    # W is the expected log return from the mistaken strategy under TRUE probabilities
    # Bet is f_B = 1/4 on B.
    # If A wins (true prob 1/2): lose bet on B. Return = 1 - 1/4 = 3/4
    # If B wins (true prob 1/4): win bet on B. Return = 1 + (1/4)*2 = 3/2
    # If C wins (true prob 1/4): lose bet on B. Return = 1 - 1/4 = 3/4
    W = p_true[0] * math.log(1 - f_mistaken_B) + \
        p_true[1] * math.log(1 + f_mistaken_B * b[1]) + \
        p_true[2] * math.log(1 - f_mistaken_B)
        
    print("Step 2: Calculate the actual growth rate W from the mistaken belief")
    print(f"The mistaken strategy is to bet f_B = {f_mistaken_B:.3f} (or {Fraction(f_mistaken_B).limit_denominator()}) on Competitor B.")
    print(f"The formula for W is: P(A)*log(Ret_A) + P(B)*log(Ret_B) + P(C)*log(Ret_C)")
    print(f"W = {p_true[0]}*log({1-f_mistaken_B}) + {p_true[1]}*log({1+f_mistaken_B*b[1]}) + {p_true[2]}*log({1-f_mistaken_B})")
    print(f"W = {W:.5f}\n")

    # Step 4: Compute the difference W* - W
    diff = W_star - W
    
    print("Step 3: Calculate the difference W* - W")
    print("The final equation for the difference is: (11/4) * log(2) - (3/2) * log(3)")
    val1 = (11/4) * math.log(2)
    val2 = (3/2) * math.log(3)
    print(f"W* - W = ({11}/{4}) * log({2}) - ({3}/{2}) * log({3})")
    print(f"W* - W = {val1:.5f} - {val2:.5f}")
    print(f"The final result is: {diff:.5f}")
    
    return diff

result = solve_growth_rate_difference()
print(f"\n<<< {result:.5f} >>>")
