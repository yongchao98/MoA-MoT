import collections

def solve_flight_revenue():
    """
    Calculates the maximum expected revenue for flight ticket sales using dynamic programming.
    """
    SEAT_CAPACITY = 10
    DAYS_TOTAL = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200

    # Demand distributions
    # {num_requests: probability}
    d1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    d2_dist = {0: 0.25, 1: 0.5, 2: 0.25}

    # dp[t][s] = max expected revenue with t days left and s seats
    dp = [[0.0 for _ in range(SEAT_CAPACITY + 1)] for _ in range(DAYS_TOTAL + 1)]
    # policy[t][s] = optimal number of cheap tickets to offer
    policy = [[0 for _ in range(SEAT_CAPACITY + 1)] for _ in range(DAYS_TOTAL + 1)]

    # Pre-compute probabilities of expensive ticket sales
    # p_sold[k][m] = {j: prob}, prob of selling j tickets to k customers with m seats available
    # max k = max d2 = 2
    p_sold = [[collections.defaultdict(float) for _ in range(SEAT_CAPACITY + 1)] for _ in range(3)]

    for m in range(SEAT_CAPACITY + 1):
        p_sold[0][m] = {0: 1.0}
    for k in range(1, 3):
        p_sold[k][0] = {0: 1.0}

    for k in range(1, 3):
        for m in range(1, SEAT_CAPACITY + 1):
            dist = collections.defaultdict(float)
            # Case: customer does not buy (prob 0.5)
            for j, p in p_sold[k-1][m].items():
                dist[j] += 0.5 * p
            # Case: customer buys (prob 0.5)
            for j, p in p_sold[k-1][m-1].items():
                dist[j + 1] += 0.5 * p
            p_sold[k][m] = dist

    # Main DP loop
    # t: days remaining
    for t in range(1, DAYS_TOTAL + 1):
        for s in range(SEAT_CAPACITY + 1):
            max_rev_for_s = -1.0
            best_c = 0
            # c: number of cheap tickets to offer
            for c in range(s + 1):
                current_c_exp_rev = 0.0
                if t <= 7:  # Week 2 (days 8-14 of sale period)
                    # Both Class 1 and Class 2 customers
                    for d1, p1 in d1_dist.items():
                        for d2, p2 in d2_dist.items():
                            sold_c2 = min(d2, c)
                            sold_c1 = min(d1, c - sold_c2)
                            
                            imm_rev_cheap = PRICE_CHEAP * (sold_c1 + sold_c2)
                            
                            d2_unhappy = d2 - sold_c2
                            s_after_cheap = s - sold_c1 - sold_c2
                            
                            dist_expensive = p_sold[d2_unhappy][s_after_cheap]
                            
                            exp_rev_expensive = 0.0
                            exp_future_rev = 0.0
                            
                            for j_sold, p_j in dist_expensive.items():
                                exp_rev_expensive += j_sold * PRICE_EXPENSIVE * p_j
                                exp_future_rev += dp[t - 1][s_after_cheap - j_sold] * p_j
                            
                            scenario_rev = imm_rev_cheap + exp_rev_expensive + exp_future_rev
                            current_c_exp_rev += p1 * p2 * scenario_rev
                else:  # Week 1 (days 1-7 of sale period)
                    # Only Class 1 customers
                    for d1, p1 in d1_dist.items():
                        sold_c = min(d1, c)
                        imm_rev = PRICE_CHEAP * sold_c
                        s_rem = s - sold_c
                        future_rev = dp[t - 1][s_rem]
                        
                        scenario_rev = imm_rev + future_rev
                        current_c_exp_rev += p1 * scenario_rev
                
                if current_c_exp_rev > max_rev_for_s:
                    max_rev_for_s = current_c_exp_rev
                    best_c = c
            
            dp[t][s] = max_rev_for_s
            policy[t][s] = best_c

    # Final result calculation and output
    max_expected_revenue = dp[DAYS_TOTAL][SEAT_CAPACITY]
    
    # Show the equation for the final step as requested
    t = DAYS_TOTAL
    s = SEAT_CAPACITY
    c_opt = policy[t][s]

    print(f"The dynamic programming solution starts {t} days before departure with {s} seats.")
    print(f"The optimal policy is to offer c = {c_opt} cheap tickets on the first day.")
    print("\nThe maximum expected revenue is calculated based on the demand distribution for that day:")
    
    # Since t=14 is in the first week, we use the simpler formula
    # E[Revenue] = E_d1[100 * sold_c + V(t-1, s - sold_c)]
    
    terms = []
    final_value = 0
    for d1, p1 in d1_dist.items():
        sold_c = min(d1, c_opt)
        imm_rev = PRICE_CHEAP * sold_c
        s_rem = s - sold_c
        future_rev = dp[t - 1][s_rem]
        term_val = imm_rev + future_rev
        terms.append({
            "prob": p1,
            "d1": d1,
            "imm_rev": imm_rev,
            "future_rev": future_rev
        })
        final_value += p1 * term_val
        
    eq_str = "Max Expected Revenue = "
    for i, term in enumerate(terms):
      eq_str += f"{term['prob']:.2f} * ({term['imm_rev']:.2f} + {term['future_rev']:.2f})"
      if i < len(terms) - 1:
        eq_str += " + "
    print(eq_str)
    
    print(f"                       = {final_value:.2f}")

    # The final answer in the requested format
    print(f"\nFinal Answer: <<<{'%.2f' % max_expected_revenue}>>>")


solve_flight_revenue()