from fractions import Fraction

def format_fraction(f):
    """Helper function to format fractions for printing."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

def main():
    """
    Solves the Kelly criterion betting problem and prints the results.
    """
    # --- Problem Setup ---
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    p_believed = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    odds = [Fraction(4), Fraction(3), Fraction(7), Fraction(7)]
    num_outcomes = len(p_true)

    # --- Step 1: Calculate the optimal bet b* and growth rate W* ---
    bets_star = [Fraction(0)] * num_outcomes
    for i in range(num_outcomes):
        # Bet only if expected payoff is greater than 1
        if p_true[i] * odds[i] > 1:
            # Kelly criterion for a single mutually exclusive outcome
            bets_star[i] = p_true[i] - (1 - p_true[i]) / (odds[i] - 1)
    
    total_bet_star = sum(bets_star)
    
    wealth_factors_star = []
    for i in range(num_outcomes):
        wealth_factors_star.append(1 - total_bet_star + bets_star[i] * odds[i])

    # --- Step 2: Calculate the incorrect bet b' based on believed probabilities ---
    bets_incorrect = [Fraction(0)] * num_outcomes
    for i in range(num_outcomes):
        if p_believed[i] * odds[i] > 1:
            bets_incorrect[i] = p_believed[i] - (1 - p_believed[i]) / (odds[i] - 1)
            
    total_bet_incorrect = sum(bets_incorrect)

    # --- Step 3: Calculate the achieved growth rate W with incorrect bet b' ---
    wealth_factors_achieved = []
    for i in range(num_outcomes):
        # Wealth factor is based on the incorrect bet
        wealth_factors_achieved.append(1 - total_bet_incorrect + bets_incorrect[i] * odds[i])

    # --- Step 4: Construct the output strings for the final equations ---

    # Full expression for W (using true probabilities and achieved wealth factors)
    w_terms = []
    for i in range(num_outcomes):
        term = f"({format_fraction(p_true[i])}) * log({format_fraction(wealth_factors_achieved[i])})"
        w_terms.append(term)
    w_expression_full = " + ".join(w_terms)

    # Simplified expression for W by grouping terms with the same wealth factor
    w_grouped_terms = {}
    for i in range(num_outcomes):
        wf = wealth_factors_achieved[i]
        w_grouped_terms[wf] = w_grouped_terms.get(wf, Fraction(0)) + p_true[i]
    w_expression_simplified = " + ".join([f"({format_fraction(coef)}) * log({format_fraction(wf)})" for wf, coef in w_grouped_terms.items()])

    # Simplified expression for W*
    # W* = (1/2)log(2) + (1/2)log(2/3) simplifies to (1/2)log(2 * 2/3) = (1/2)log(4/3)
    w_star_expression_simplified = f"({format_fraction(Fraction(1, 2))}) * log({format_fraction(Fraction(4, 3))})"

    # --- Print the results ---
    print("If you bet based on your incorrect probabilities, the doubling rate W you will achieve is:")
    print(f"W = {w_expression_full}")
    print("\nSimplifying by combining terms with the same outcome, the equation is:")
    print(f"W = {w_expression_simplified}")
    print("\nThe decrease in your doubling rate, Delta W = W* - W, is:")
    delta_w_expression = f"{w_star_expression_simplified} - ({w_expression_simplified})"
    print(f"ΔW = {delta_w_expression}")

if __name__ == "__main__":
    main()