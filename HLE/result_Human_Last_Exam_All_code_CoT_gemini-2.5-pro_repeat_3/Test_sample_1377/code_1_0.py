import math

def solve_race_growth_rate():
    """
    Solves the Kelly Criterion problem as described.
    Calculates W* (optimal growth rate) and W (actual growth rate from mistaken beliefs)
    and returns their difference.
    """
    # Use natural logarithm as is standard for Kelly Criterion
    log = math.log

    # -- Step 1: Define Inputs --
    # Actual probabilities
    p_actual = [1/2, 1/4, 1/4]
    # Believed probabilities
    p_believed = [1/4, 1/2, 1/4]
    # Net odds (payout ratio b+1:1 means net odds are b)
    # 4:1 payout -> b=3. 3:1 payout -> b=2
    b = [3, 2, 2]
    competitors = ['A', 'B', 'C']

    # -- Step 2: Calculate Optimal Strategy (f*) --
    f_star = [0, 0, 0]
    # Check expectation E = p*b - (1-p) to see if a bet is favorable
    E_actual = [p_actual[i] * b[i] - (1 - p_actual[i]) for i in range(3)]
    
    # Simple Kelly fraction for single positive expectation bet
    # f = E/b. We bet only if E > 0.
    for i in range(3):
        if E_actual[i] > 0:
            f_star[i] = E_actual[i] / b[i]

    # -- Step 3: Calculate Optimal Growth Rate (W*) --
    # Wealth outcomes: W_i = 1 - sum(f) + f_i*(b_i+1) = 1 + f_i*b_i - sum(f_j for j!=i)
    total_f_star = sum(f_star)
    wealth_outcomes_star = [1 - total_f_star + f_star[i] * (b[i] + 1) for i in range(3)]
    
    W_star = sum(p_actual[i] * log(wealth_outcomes_star[i]) for i in range(3))

    # -- Step 4: Calculate Mistaken Strategy (f') --
    f_mistaken = [0, 0, 0]
    # Check expectation E' = p'*b - (1-p') with believed probabilities
    E_believed = [p_believed[i] * b[i] - (1 - p_believed[i]) for i in range(3)]

    # Calculate mistaken fractions based on perceived positive expectation
    for i in range(3):
        if E_believed[i] > 0:
            f_mistaken[i] = E_believed[i] / b[i]

    # -- Step 5: Calculate Actual Growth Rate (W) --
    # Calculate wealth outcomes using mistaken fractions
    total_f_mistaken = sum(f_mistaken)
    wealth_outcomes_mistaken = [1 - total_f_mistaken + f_mistaken[i] * (b[i] + 1) for i in range(3)]

    # Calculate actual growth rate W using actual probabilities and mistaken fractions
    W_actual = sum(p_actual[i] * log(wealth_outcomes_mistaken[i]) for i in range(3))
    
    # -- Step 6: Calculate the Difference and Print --
    difference = W_star - W_actual
    
    # The analytical result is (11/4)*log(2) - (3/2)*log(3)
    c1 = 11/4
    v1 = 2
    c2 = 3/2
    v2 = 3
    
    print("The final result W* - W is calculated from the equation:")
    print(f"({c1}) * log({v1}) - ({c2}) * log({v2})")
    
    print("\nSubstituting the values:")
    term1_val = c1 * log(v1)
    term2_val = c2 * log(v2)
    print(f"= ({c1}) * {log(v1):.5f} - ({c2}) * {log(v2):.5f}")
    print(f"= {term1_val:.5f} - {term2_val:.5f}")
    
    print("\nThe final numerical value for W* - W is:")
    print(f"{difference:.8f}")

solve_race_growth_rate()
<<<0.25823621>>>