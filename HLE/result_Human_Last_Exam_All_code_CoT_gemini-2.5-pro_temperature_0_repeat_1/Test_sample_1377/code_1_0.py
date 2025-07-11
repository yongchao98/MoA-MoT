import math

def calculate_growth_difference():
    """
    Calculates the difference between the optimal and actual wealth growth rates
    in a betting scenario with mistaken probabilities.
    """
    # --- Problem Parameters ---
    # True probabilities
    p = {'A': 1/2, 'B': 1/4, 'C': 1/4}
    # Mistaken probabilities
    q = {'A': 1/4, 'B': 1/2, 'C': 1/4}
    # Payout ratios (b:1), which means net odds are b
    b = {'A': 3, 'B': 2, 'C': 2}

    print("This script calculates the difference between optimal (W*) and actual (W) growth rates.\n")

    # --- 1. Calculate Optimal Growth Rate (W*) ---
    print("--- Step 1: Calculating the Optimal Growth Rate (W*) ---")
    # The Kelly criterion suggests betting only when p*(b+1) > 1.
    # For the true probabilities:
    # A: (1/2) * (3 + 1) = 2.0 > 1 (Bet)
    # B: (1/4) * (2 + 1) = 0.75 < 1 (Don't bet)
    # C: (1/4) * (2 + 1) = 0.75 < 1 (Don't bet)
    # So, the optimal strategy is to only bet on A.

    # Optimal fraction to bet on A using the Kelly formula: f = (p*(b+1) - 1) / b
    f_A_star = (p['A'] * (b['A'] + 1) - 1) / b['A']

    # Optimal growth rate W* = p_A*log(1+f_A*b_A) + (p_B+p_C)*log(1-f_A)
    w_star = p['A'] * math.log(1 + f_A_star * b['A']) + \
             (p['B'] + p['C']) * math.log(1 - f_A_star)

    print("The optimal strategy is to bet a fraction of 1/3 on Competitor A.")
    print("The equation for W* is: p(A)*log(1 + f_A*b_A) + (p(B)+p(C))*log(1 - f_A)")
    print(f"W* = ({p['A']}) * log(1 + (1/3) * {b['A']}) + ({p['B']} + {p['C']}) * log(1 - (1/3))")
    print(f"W* = {p['A']} * log({1 + f_A_star * b['A']:.4f}) + {p['B'] + p['C']} * log({1 - f_A_star:.4f})")
    print(f"W* = {w_star:.4f}\n")

    # --- 2. Calculate Actual Growth Rate (W) from Mistaken Strategy ---
    print("--- Step 2: Calculating the Actual Growth Rate (W) ---")
    # Based on mistaken probabilities q:
    # A: (1/4) * (3 + 1) = 1.0 (Don't bet, edge is zero)
    # B: (1/2) * (2 + 1) = 1.5 > 1 (Bet)
    # C: (1/4) * (2 + 1) = 0.75 < 1 (Don't bet)
    # So, the mistaken strategy is to only bet on B.

    # Mistaken fraction to bet on B: f = (q*(b+1) - 1) / b
    f_B_mistake = (q['B'] * (b['B'] + 1) - 1) / b['B']

    # Actual growth rate W is calculated using f_B_mistake but with true probabilities p.
    # W = p_A*log(1-f_B) + p_B*log(1+f_B*b_B) + p_C*log(1-f_B)
    w = (p['A'] + p['C']) * math.log(1 - f_B_mistake) + \
        p['B'] * math.log(1 + f_B_mistake * b['B'])

    print("The mistaken strategy is to bet a fraction of 1/4 on Competitor B.")
    print("The equation for W is: (p(A)+p(C))*log(1 - f_B) + p(B)*log(1 + f_B*b_B)")
    print(f"W = ({p['A']} + {p['C']}) * log(1 - (1/4)) + ({p['B']}) * log(1 + (1/4) * {b['B']})")
    print(f"W = {p['A'] + p['C']} * log({1 - f_B_mistake:.4f}) + {p['B']} * log({1 + f_B_mistake * b['B']:.4f})")
    print(f"W = {w:.4f}\n")

    # --- 3. Calculate the Difference W* - W ---
    print("--- Step 3: Calculating the Final Difference (W* - W) ---")
    difference = w_star - w
    print("The final equation is W* - W:")
    print(f"{w_star:.4f} - ({w:.4f}) = {difference:.4f}")

    return difference

# Run the calculation and store the final answer
final_answer = calculate_growth_difference()
# The final answer is returned in the specified format below.
# print(f"\n<<< {final_answer:.4f} >>>")