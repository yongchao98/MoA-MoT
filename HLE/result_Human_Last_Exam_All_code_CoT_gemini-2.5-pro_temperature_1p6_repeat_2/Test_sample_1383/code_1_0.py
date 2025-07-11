import math
from fractions import Fraction

def solve_betting_problem():
    """
    Calculates and explains the solution to the Kelly criterion betting problem.
    """

    # Helper function to format fractions for clean printing
    def format_frac(f):
        if f.denominator == 1:
            return str(f.numerator)
        return f"{f.numerator}/{f.denominator}"

    # --- Problem Setup ---
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    p_false = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    o = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]  # Total return odds (X-for-1)
    d = [x - 1 for x in o]  # Net payoff odds
    num_bikes = len(p_true)
    
    # --- PART 1: Find W (Achieved Growth Rate with incorrect probabilities) ---

    print("Step 1: Determine the betting strategy based on your incorrect probabilities (q).")
    print("We use the Kelly criterion, which states a bet should be made only if the expected value is positive (q_i * o_i > 1).\n")
    
    q = p_false
    b_false = [Fraction(0)] * num_bikes
    bet_indices_false = []
    
    for i in range(num_bikes):
        val = q[i] * o[i]
        decision = "> 1, so a bet is warranted." if val > 1 else "<= 1, so no bet should be placed."
        print(f"Bike {i+1}: q_{i+1}*o_{i+1} = {format_frac(q[i])} * {format_frac(o[i])} = {format_frac(val)}. This is {decision}")
        if val > 1:
            bet_indices_false.append(i)

    # In this case, only Bike 2 meets the criterion.
    i = bet_indices_false[0]
    q_i = q[i]
    d_i = d[i]
    b_false[i] = q_i - (1 - q_i) / d_i
    
    print(f"\nBased on your beliefs, you only bet on Bike {i+1}. The Kelly fraction is b = q - (1-q)/d.")
    print(f"b_2 = q_2 - (1 - q_2) / d_2")
    print(f"b_2 = {format_frac(q_i)} - (1 - {format_frac(q_i)}) / {format_frac(d_i)}")
    print(f"b_2 = {format_frac(q_i)} - {format_frac(1-q_i)} / {format_frac(d_i)} = {format_frac(q_i)} - {format_frac((1-q_i)/d_i)} = {format_frac(b_false[i])}")
    print(f"So, your strategy is to bet {format_frac(b_false[i])} of your wealth on Bike 2.")

    print("\nStep 2: Calculate the actually achieved growth rate (W) using this strategy and the TRUE probabilities (p).")
    print("W = Σ p_i * ln(Wealth Factor if bike i wins)")
    
    total_bet_false = sum(b_false)
    w_terms = [{'p': p_true[i], 'growth': (1 - total_bet_false) + b_false[i] * o[i]} for i in range(num_bikes)]

    # Group terms for the final equation of W
    grouped_terms_w = {}
    for term in w_terms:
        coeff, log_arg = term['p'], term['growth']
        if log_arg in grouped_terms_w:
            grouped_terms_w[log_arg] += coeff
        else:
            grouped_terms_w[log_arg] = coeff
            
    w_str_parts = []
    # Sort for consistent output
    for log_arg, coeff in sorted(grouped_terms_w.items(), key=lambda item: item[0]):
        w_str_parts.append(f"({format_frac(coeff)}) * ln({format_frac(log_arg)})")
    
    print("\nThe achieved doubling rate W is:")
    final_W_expression = " + ".join(w_str_parts)
    print("W = " + final_W_expression)
    
    # --- PART 2: Find ΔW = W* - W ---
    
    print("\n" + "="*50)
    print("\nStep 3: Calculate the optimal growth rate W* (using true probabilities for the strategy).")
    
    # Determine optimal bet b* based on true probabilities p_true
    b_star = [Fraction(0)] * num_bikes
    i_star = [i for i in range(num_bikes) if p_true[i] * o[i] > 1][0]
    b_star[i_star] = p_true[i_star] - (1 - p_true[i_star]) / d[i_star]
    
    # Calculate W*
    total_bet_star = sum(b_star)
    w_star_terms = [{'p': p_true[i], 'growth': (1 - total_bet_star) + b_star[i] * o[i]} for i in range(num_bikes)]
        
    grouped_terms_star = {}
    for term in w_star_terms:
        coeff, log_arg = term['p'], term['growth']
        if log_arg in grouped_terms_star: grouped_terms_star[log_arg] += coeff
        else: grouped_terms_star[log_arg] = coeff
            
    w_star_str_parts = []
    for log_arg, coeff in sorted(grouped_terms_star.items(), key=lambda item: item[0]):
        w_star_str_parts.append(f"({format_frac(coeff)}) * ln({format_frac(log_arg)})")
    
    W_star_final_expression = " + ".join(w_star_str_parts)
    print("Following the same logic, the truly optimal strategy is to bet 1/3 on Bike 1.")
    print(f"The resulting optimal growth rate W* is: W* = {W_star_final_expression}")

    print("\nStep 4: Calculate the decrease in growth rate, ΔW = W* - W.")
    print("To simplify the subtraction, we can expand the logarithms: ln(a/b) = ln(a) - ln(b).")
    # Manually derived simplified forms
    print("W* = ln(2) - (1/2) * ln(3)")
    print("W  = (3/4) * ln(3/4) + (1/4) * ln(3/2) = ln(3) - (7/4) * ln(2)")
    print("\nΔW = W* - W = (ln(2) - (1/2)ln(3)) - (ln(3) - (7/4)ln(2))")
    
    # Final coefficients for ΔW = c1*ln(2) - c2*ln(3)
    c1 = Fraction(11, 4)
    c2 = Fraction(3, 2)
    print(f"\nTherefore, the decrease in the doubling rate is:")
    print(f"ΔW = ({format_frac(c1)}) * ln(2) - ({format_frac(c2)}) * ln(3)")

if __name__ == '__main__':
    solve_betting_problem()
