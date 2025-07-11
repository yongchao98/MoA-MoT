import collections

def solve_flight_revenue():
    """
    Calculates the maximum expected total revenue for the flight using dynamic programming.
    """
    # --- Problem Constants ---
    CAPACITY = 10
    DAYS = 14
    CHEAP_PRICE = 100
    EXPENSIVE_PRICE = 200

    # --- Customer Demand Distributions ---
    DIST_C1 = {0: 0.25, 1: 0.5, 2: 0.25}
    DIST_C2_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25} # Days 8-14 of selling period (t<=7)
    DIST_C2_WEEK1 = {0: 1.0, 1: 0.0, 2: 0.0}   # Days 1-7 of selling period (t>7)

    # --- Binomial Probabilities P(X=k) for X~B(n, 0.5), for n up to 2 ---
    BINOM_PROBS = {
        0: {0: 1.0},
        1: {0: 0.5, 1: 0.5},
        2: {0: 0.25, 1: 0.5, 2: 0.25},
    }

    def calculate_expected_revenue_for_decision(t, s, c, v_prev):
        """
        Calculates the expected total revenue for a state (t, s) and decision c.
        t: days remaining, s: seats available, c: cheap tickets to offer
        v_prev: DP table for t-1 days remaining (V[t-1])
        """
        dist_c2 = DIST_C2_WEEK2 if t <= 7 else DIST_C2_WEEK1
        total_expected_revenue = 0
        
        for n1, p1 in DIST_C1.items():
            for n2, p2 in dist_c2.items():
                
                # --- Cheap ticket sales (Class 2 has priority) ---
                sold_c2 = min(n2, c)
                sold_c1 = min(n1, c - sold_c2)
                sold_cheap_total = sold_c1 + sold_c2
                revenue_from_cheap = sold_cheap_total * CHEAP_PRICE
                
                # --- Expensive ticket sales (probabilistic) ---
                unfulfilled_c2 = n2 - sold_c2
                expensive_seats_available = s - c
                
                exp_value_over_k = 0
                if unfulfilled_c2 == 0 or expensive_seats_available == 0:
                    seats_left = s - sold_cheap_total
                    exp_value_over_k = revenue_from_cheap + v_prev[seats_left]
                else:
                    for k, prob_k in BINOM_PROBS[unfulfilled_c2].items():
                        sold_expensive = min(k, expensive_seats_available)
                        revenue_from_expensive = sold_expensive * EXPENSIVE_PRICE
                        seats_left = s - sold_cheap_total - sold_expensive
                        
                        total_rev_outcome = revenue_from_cheap + revenue_from_expensive + v_prev[seats_left]
                        exp_value_over_k += prob_k * total_rev_outcome
                        
                total_expected_revenue += p1 * p2 * exp_value_over_k
                
        return total_expected_revenue

    # --- Main DP Calculation ---
    V = collections.defaultdict(lambda: collections.defaultdict(float))

    for t in range(1, DAYS + 1):
        for s in range(CAPACITY + 1):
            revenues_for_each_c = [calculate_expected_revenue_for_decision(t, s, c, V[t-1]) for c in range(s + 1)]
            V[t][s] = max(revenues_for_each_c) if revenues_for_each_c else 0

    # --- Output the Explanation and Final Result ---
    print("This problem is solved using dynamic programming by calculating the maximum expected revenue V(t, s) for each state (t=days left, s=seats left).")
    print("The final result V(14, 10) is found by selecting the best decision 'c' (number of cheap tickets to offer) on the first day.")
    print("\nBelow is the breakdown of the expected total revenue for each possible decision 'c' at the start (t=14, s=10):")
    
    final_revenues = []
    for c in range(CAPACITY + 1):
        rev = calculate_expected_revenue_for_decision(DAYS, CAPACITY, c, V[DAYS-1])
        final_revenues.append(rev)
        # The equation for each E[c] is an extensive sum over all possibilities.
        # Here we just show the resulting value for each top-level decision.
        print(f"  Decision c = {c:2d} (offer {c} cheap tickets): Expected Revenue = {rev:.2f}")

    max_revenue = V[DAYS][CAPACITY]
    optimal_c = final_revenues.index(max_revenue)

    print(f"\nThe optimal decision on the first day is c = {optimal_c}, which yields the maximum expected revenue.")
    
    # Building the "final equation" string as requested
    revenue_strings = [f"{rev:.2f}" for rev in final_revenues]
    final_equation = f"Max Expected Revenue = max({', '.join(revenue_strings)})"
    
    print("\nThe final equation is:")
    print(f"{final_equation} = {max_revenue:.4f}")

    print(f"\nMaximum expected total revenue:")
    print(f"<<<{max_revenue:.4f}>>>")

solve_flight_revenue()