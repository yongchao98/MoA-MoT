def solve_bike_race_problem():
    """
    This script calculates the achieved growth rate (W) and the decrease in
    growth rate (Delta W) for the bike race betting problem. It presents the
    answers in terms of natural logs and fractions.
    """

    # --- Problem Data ---
    # True probabilities
    p = [1/2, 1/4, 1/8, 1/8]
    p_str = ["1/2", "1/4", "1/8", "1/8"]
    # Incorrectly believed probabilities (used for betting fractions)
    q = [1/4, 1/2, 1/8, 1/8]
    q_str = ["1/4", "1/2", "1/8", "1/8"]
    # Bookmaker's odds (o-for-1)
    o = [4, 3, 7, 7]
    o_str = ["4", "3", "7", "7"]

    # --- Part 1: Calculate the achieved growth rate W ---
    print("--- Achieved Growth Rate (W) with incorrect beliefs ---")
    print("The formula is W = sum(p_i * log(q_i * o_i)) where q_i are the betting fractions.")
    
    # Build the full equation string
    w_eq_full = "W = "
    w_eq_terms = []
    for i in range(len(p)):
        w_eq_terms.append(f"({p_str[i]}) * log(({q_str[i]}) * {o_str[i]})")
    w_eq_full += " + ".join(w_eq_terms)
    print(w_eq_full)

    # Build the simplified equation string
    w_eq_simplified = "W = "
    w_eq_simp_terms = [
        "(1/2) * log(1)",
        "(1/4) * log(3/2)",
        "(1/8) * log(7/8)",
        "(1/8) * log(7/8)"
    ]
    w_eq_simplified += " + ".join(w_eq_simp_terms)
    print(w_eq_simplified)
    
    # Final expression for W
    w_final = "W = (1/4) * log(3/2) + (1/4) * log(7/8)"
    print("Since log(1) = 0, this simplifies to:")
    print(w_final)
    print("-" * 20)

    # --- Part 2: Calculate the decrease in growth rate Delta W ---
    print("\n--- Decrease in Growth Rate (Delta W) ---")
    print("To find the decrease, we first need the optimal rate W*.")
    
    # Calculation of W*
    print("W* is calculated by betting the true probabilities: W* = sum(p_i * log(p_i * o_i)).")
    w_star_eq = "W* = (1/2) * log((1/2) * 4) + (1/4) * log((1/4) * 3) + (1/8) * log((1/8) * 7) + (1/8) * log((1/8) * 7)"
    print(w_star_eq)
    w_star_final = "W* = (1/2) * log(2) + (1/4) * log(3/4) + (1/4) * log(7/8)"
    print("Simplifying gives:")
    print(w_star_final)
    
    # Calculation of Delta W = W* - W
    print("\nNow we calculate Delta W = W* - W:")
    delta_w_eq = f"Delta W = ({w_star_final.replace('W* = ','')}) - ({w_final.replace('W = ','')})"
    print(delta_w_eq)

    # Simplifying Delta W
    delta_w_simplified = "Delta W = (1/2) * log(2) + (1/4) * log(3/4) - (1/4) * log(3/2)"
    print("The term (1/4) * log(7/8) cancels out, leaving:")
    print(delta_w_simplified)

    # Final expression for Delta W
    delta_w_final = "Delta W = (1/4) * log(2)"
    print("Using log properties (log(a/b) = log(a)-log(b)), this simplifies to:")
    print(delta_w_final)
    print("-" * 20)

if __name__ == '__main__':
    solve_bike_race_problem()