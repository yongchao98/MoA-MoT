import sympy

def solve_bike_race_betting():
    """
    This function calculates the optimal and achieved logarithmic growth rates
    for a betting scenario and the difference between them.
    """
    # Define probabilities and odds using SymPy's Rational for precision
    p_true = [sympy.Rational(1, 2), sympy.Rational(1, 4), sympy.Rational(1, 8), sympy.Rational(1, 8)]
    q_believed = [sympy.Rational(1, 4), sympy.Rational(1, 2), sympy.Rational(1, 8), sympy.Rational(1, 8)]
    # "x-for-1" odds correspond to decimal odds of x
    odds = [sympy.Rational(o) for o in [4, 3, 7, 7]]

    # --- 1. Calculate the Optimal Growth Rate W* ---
    # The optimal strategy is to bet fractions b_i = p_i.
    # The formula is W* = sum(p_i * ln(p_i * o_i)).
    W_star = sum(p * sympy.log(p * o) for p, o in zip(p_true, odds))

    print("1. Calculating the optimal growth rate W*:")
    print("W* = sum(p_i * ln(p_i * o_i))")
    
    # Show the full equation with numbers
    w_star_eq_parts = []
    for p, o in zip(p_true, odds):
        w_star_eq_parts.append(f"({p}) * ln(({p}) * ({o}))")
    print(f"W* = {' + '.join(w_star_eq_parts)}")

    # Show the equation with simplified terms inside the log
    w_star_simplified_terms = []
    for p, o in zip(p_true, odds):
         w_star_simplified_terms.append(f"({p}) * ln({p * o})")
    print(f"W* = {' + '.join(w_star_simplified_terms)}")

    # Calculate and print the final symbolic expression for W*
    W_star_final = sympy.expand_log(W_star, force=True).simplify()
    print(f"\nSimplified W* = {W_star_final}")
    print("-" * 50)

    # --- 2. Calculate the Achieved Growth Rate W ---
    # This is calculated by betting fractions b_i = q_i, while true probabilities are p_i.
    # The formula is W = sum(p_i * ln(q_i * o_i)).
    W = sum(p * sympy.log(q * o) for p, q, o in zip(p_true, q_believed, odds))
    
    print("2. Calculating the achieved growth rate W (based on incorrect beliefs):")
    print("W = sum(p_i * ln(q_i * o_i))")
    
    # Show the full equation with numbers
    w_eq_parts = []
    for p, q, o in zip(p_true, q_believed, odds):
        w_eq_parts.append(f"({p}) * ln(({q}) * ({o}))")
    print(f"W = {' + '.join(w_eq_parts)}")

    # Show the equation with simplified terms inside the log
    w_simplified_terms = []
    for p, q, o in zip(p_true, q_believed, odds):
        w_simplified_terms.append(f"({p}) * ln({q * o})")
    print(f"W = {' + '.join(w_simplified_terms)}")
    
    # Calculate and print the final symbolic expression for W
    W_final = sympy.expand_log(W, force=True).simplify()
    print(f"\nThe achieved doubling rate is W = {W_final}")
    print("-" * 50)

    # --- 3. Calculate the Decrease in Growth Rate Delta W ---
    # Delta W = W* - W, which is also the KL divergence D_KL(p || q).
    # Delta W = sum(p_i * ln(p_i / q_i))
    Delta_W = sum(p * sympy.log(p / q) for p, q in zip(p_true, q_believed))
    
    print("3. Calculating the decrease in growth rate Delta W = W* - W:")
    print("This is equivalent to the Kullback-Leibler divergence D_KL(p || q).")
    print("Delta W = sum(p_i * ln(p_i / q_i))")
    
    # Show the full equation with numbers
    delta_w_eq_parts = []
    for p, q in zip(p_true, q_believed):
        delta_w_eq_parts.append(f"({p}) * ln(({p}) / ({q}))")
    print(f"Delta W = {' + '.join(delta_w_eq_parts)}")
    
    # Show the equation with simplified terms inside the log
    delta_w_simplified_terms = []
    for p, q in zip(p_true, q_believed):
        delta_w_simplified_terms.append(f"({p}) * ln({p / q})")
    print(f"Delta W = {' + '.join(delta_w_simplified_terms)}")

    # Calculate and print the final symbolic expression for Delta W
    Delta_W_final = sympy.expand_log(Delta_W, force=True).simplify()
    print(f"\nThe decrease in doubling rate is Delta W = {Delta_W_final}")
    print("-" * 50)

if __name__ == '__main__':
    solve_bike_race_betting()