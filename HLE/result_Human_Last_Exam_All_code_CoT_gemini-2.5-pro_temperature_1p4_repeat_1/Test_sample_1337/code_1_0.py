import math

# Use a memoization dictionary for the combinations function to speed it up.
comb_memo = {}

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    """
    if k < 0 or k > n:
        return 0
    if (n, k) in comb_memo:
        return comb_memo[(n, k)]
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    comb_memo[(n, k)] = res
    return res

def solve():
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    # --- Problem Parameters ---
    TOTAL_SEATS = 10
    NUM_DAYS = 14
    CHEAP_PRICE = 100
    EXPENSIVE_PRICE = 200

    # Customer demand distributions
    P1 = {0: 0.25, 1: 0.5, 2: 0.25}
    P2_active = {0: 0.25, 1: 0.5, 2: 0.25}
    P2_inactive = {0: 1.0}

    # --- DP Initialization ---
    max_expected_revenue = -1.0
    best_protection_level = -1
    best_E_table = None

    # --- Iterate through all possible protection levels (policies) ---
    for p in range(TOTAL_SEATS + 1):
        c_limit = TOTAL_SEATS - p
        e_limit = p

        # DP table: E[day][cheap_seats][expensive_seats]
        E = [[[0.0 for _ in range(e_limit + 1)] for _ in range(c_limit + 1)] for _ in range(NUM_DAYS + 2)]
        
        # --- DP Calculation (backward induction) ---
        for t in range(NUM_DAYS, 0, -1):
            P2 = P2_active if t > 7 else P2_inactive
            
            for c in range(c_limit + 1):
                for e in range(e_limit + 1):
                    
                    total_expected_value = 0.0
                    
                    # Iterate over all possible demand scenarios for day t
                    for d1, prob1 in P1.items():
                        for d2, prob2 in P2.items():
                            prob_scenario = prob1 * prob2
                            if prob_scenario == 0:
                                continue

                            # -- Simulate sales for the scenario (d1, d2) --
                            
                            # 1. Class 2 customers purchase cheap tickets (they have priority)
                            cheap_sales_c2 = min(c, d2)
                            rev_cheap_c2 = cheap_sales_c2 * CHEAP_PRICE
                            c_after_c2 = c - cheap_sales_c2
                            d2_unfulfilled = d2 - cheap_sales_c2

                            # 2. Class 1 customers purchase cheap tickets
                            cheap_sales_c1 = min(c_after_c2, d1)
                            rev_cheap_c1 = cheap_sales_c1 * CHEAP_PRICE
                            c_final = c_after_c2 - cheap_sales_c1

                            # 3. Class 2 customers consider expensive tickets
                            exp_rev_exp_and_future = 0.0
                            if d2_unfulfilled > 0 and e > 0:
                                prob_buy_exp = 0.5
                                # Iterate over k, the number of C2 customers who try to buy an expensive ticket
                                for k in range(d2_unfulfilled + 1):
                                    prob_k_buy = combinations(d2_unfulfilled, k) * (prob_buy_exp**k) * ((1 - prob_buy_exp)**(d2_unfulfilled - k))
                                    
                                    expensive_sales = min(e, k)
                                    rev_exp = expensive_sales * EXPENSIVE_PRICE
                                    e_final = e - expensive_sales
                                    
                                    future_rev = E[t + 1][c_final][e_final]
                                    exp_rev_exp_and_future += prob_k_buy * (rev_exp + future_rev)
                            else:
                                # No unfulfilled C2 customers or no expensive seats left for them
                                e_final = e
                                future_rev = E[t + 1][c_final][e_final]
                                exp_rev_exp_and_future = future_rev
                            
                            scenario_value = rev_cheap_c1 + rev_cheap_c2 + exp_rev_exp_and_future
                            total_expected_value += prob_scenario * scenario_value
                    
                    E[t][c][e] = total_expected_value

        current_revenue = E[1][c_limit][e_limit]
        if current_revenue > max_expected_revenue:
            max_expected_revenue = current_revenue
            best_protection_level = p
            best_E_table = E

    # --- Construct and Print Final Equation ---
    p = best_protection_level
    c_initial = TOTAL_SEATS - p
    e_initial = p
    
    # Calculate components for the Day 1 expectation equation
    # On Day 1, only Class 1 customers arrive (d2=0)
    
    # Case d1 = 0
    prob_d1_0 = P1[0]
    rev_d1_0 = min(c_initial, 0) * CHEAP_PRICE
    c_next_0 = c_initial - min(c_initial, 0)
    future_rev_0 = best_E_table[2][c_next_0][e_initial]
    
    # Case d1 = 1
    prob_d1_1 = P1[1]
    rev_d1_1 = min(c_initial, 1) * CHEAP_PRICE
    c_next_1 = c_initial - min(c_initial, 1)
    future_rev_1 = best_E_table[2][c_next_1][e_initial]

    # Case d1 = 2
    prob_d1_2 = P1[2]
    rev_d1_2 = min(c_initial, 2) * CHEAP_PRICE
    c_next_2 = c_initial - min(c_initial, 2)
    future_rev_2 = best_E_table[2][c_next_2][e_initial]

    # Print the equation for the first day's expected revenue calculation
    # which results in the maximum total expected revenue
    print(
        f"{prob_d1_0:.2f} * ({rev_d1_0} + {future_rev_0:.2f}) + "
        f"{prob_d1_1:.2f} * ({rev_d1_1} + {future_rev_1:.2f}) + "
        f"{prob_d1_2:.2f} * ({rev_d1_2} + {future_rev_2:.2f}) = {max_expected_revenue:.2f}"
    )

solve()
<<<1276.54>>>