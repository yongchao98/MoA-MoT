import math

def solve_growth_rates():
    """
    Calculates the optimal vs actual growth rate for the given race scenario.
    """
    # 1. Define probabilities and payouts
    # True probabilities
    p_true = {'A': 0.5, 'B': 0.25, 'C': 0.25}
    # Mistaken probabilities
    p_mistaken = {'A': 0.25, 'B': 0.5, 'C': 0.25}
    # Payout multipliers (odds + 1)
    odds = {'A': 5, 'B': 4, 'C': 4}

    # 2. Calculate Optimal Fractions and W* (based on true probabilities)
    # The optimal fractions f*_k are found by solving p_k*o_k / R_k = constant
    # which leads to f*_A=0.5, f*_B=0.25, f*_C=0.25
    f_star = {'A': 0.5, 'B': 0.25, 'C': 0.25}

    # Calculate W*
    # W* = p_A*ln(o_A*f*_A) + p_B*ln(o_B*f*_B) + p_C*ln(o_C*f*_C)
    # Here total fraction is 1. R_k simplifies to o_k * f_k
    w_star = (p_true['A'] * math.log(odds['A'] * f_star['A']) +
              p_true['B'] * math.log(odds['B'] * f_star['B']) +
              p_true['C'] * math.log(odds['C'] * f_star['C']))
    
    # 3. Calculate Mistaken Fractions and W (based on mistaken probabilities)
    # The mistaken fractions f'_k are found by solving p'_k*o_k / R'_k = constant
    # which leads to f'_A=0.25, f'_B=0.5, f'_C=0.25
    f_mistaken = {'A': 0.25, 'B': 0.5, 'C': 0.25}

    # Calculate W using mistaken fractions but true probabilities
    # W = p_A*ln(o_A*f'_A) + p_B*ln(o_B*f'_B) + p_C*ln(o_C*f'_C)
    w_actual = (p_true['A'] * math.log(odds['A'] * f_mistaken['A']) +
                p_true['B'] * math.log(odds['B'] * f_mistaken['B']) +
                p_true['C'] * math.log(odds['C'] * f_mistaken['C']))
    
    # 4. Calculate the difference
    difference = w_star - w_actual

    # Print the equation with calculated numbers
    print(f"Optimal Growth Rate (W*): {w_star:.8f}")
    print(f"Actual Growth Rate (W): {w_actual:.8f}")
    print(f"Difference (W* - W): {w_star:.8f} - {w_actual:.8f} = {difference:.8f}")

    # A more direct calculation from the analytical solution
    # difference = 0.25 * ln(2)
    analytical_diff = 0.25 * math.log(2)
    print(f"\nAnalytical solution: 0.25 * ln(2) = {analytical_diff:.8f}")

solve_growth_rates()
