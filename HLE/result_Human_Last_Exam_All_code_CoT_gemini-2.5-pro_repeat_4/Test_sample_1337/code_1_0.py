import functools

def solve_flight_revenue():
    """
    This function calculates the maximum expected revenue for selling flight tickets
    by finding the optimal protection level 'k'.
    """
    
    # --- Problem Parameters ---
    SEAT_CAPACITY = 10
    NUM_DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200
    PROB_BUY_EXPENSIVE = 0.5
    
    # Customer arrival distributions
    # P(n requests) = dist[n]
    P1_DIST = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1 customers (every day)
    P2_DIST = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2 customers (second week only)

    max_revenue = -1.0
    best_k = -1
    all_k_revenues = {}

    # --- Part 1: Find the optimal protection level k ---
    # A protection level 'k' means cheap tickets are only sold if more than k seats are available.
    
    print("Step 1: Calculating expected revenue for each possible protection level (k)...")

    for k in range(SEAT_CAPACITY + 1):
        # DP table: E[t][s] = Expected revenue with t days remaining and s seats available.
        E = [[0.0] * (SEAT_CAPACITY + 1) for _ in range(NUM_DAYS + 1)]

        for t in range(1, NUM_DAYS + 1):  # t = days remaining until departure
            day_index = NUM_DAYS - t + 1  # Day number from 1 to 14
            is_second_week = day_index > 7

            # This recursive function calculates the total expected value for a given day's arrivals.
            # It's defined inside the loop to capture the correct 't' and 'k'.
            @functools.lru_cache(maxsize=None)
            def get_daily_total_value(n1, n2, s_avail):
                if s_avail <= 0 or (n1 == 0 and n2 == 0):
                    return E[t - 1][max(0, s_avail)]

                if n2 > 0:  # Class 2 has priority
                    if s_avail > k:
                        # Sell cheap ticket ($100)
                        return PRICE_CHEAP + get_daily_total_value(n1, n2 - 1, s_avail - 1)
                    else:
                        # Offer expensive ticket ($200)
                        val_if_buy = PRICE_EXPENSIVE + get_daily_total_value(n1, n2 - 1, s_avail - 1)
                        val_if_walk = get_daily_total_value(n1, n2 - 1, s_avail)
                        return PROB_BUY_EXPENSIVE * val_if_buy + (1 - PROB_BUY_EXPENSIVE) * val_if_walk
                
                if n1 > 0:  # Process Class 1
                    if s_avail > k:
                        # Sell cheap ticket ($100)
                        return PRICE_CHEAP + get_daily_total_value(n1 - 1, n2, s_avail - 1)
                    else:
                        # Class 1 walks away
                        return get_daily_total_value(n1 - 1, n2, s_avail)
                return E[t - 1][s_avail] # Should not be reached if n1/n2>0

            for s in range(1, SEAT_CAPACITY + 1):
                exp_rev_for_ts = 0.0
                p2_dist_today = P2_DIST if is_second_week else {0: 1.0}

                for n1_arrivals, p1 in P1_DIST.items():
                    for n2_arrivals, p2 in p2_dist_today.items():
                        prob_arrival = p1 * p2
                        total_value = get_daily_total_value(n1_arrivals, n2_arrivals, s)
                        exp_rev_for_ts += prob_arrival * total_value
                
                E[t][s] = exp_rev_for_ts
            
            get_daily_total_value.cache_clear() # Clear cache for the next day 't'

        final_revenue_for_k = E[NUM_DAYS][SEAT_CAPACITY]
        all_k_revenues[k] = final_revenue_for_k

        if final_revenue_for_k > max_revenue:
            max_revenue = final_revenue_for_k
            best_k = k
            
    print("Calculation complete.\n")
    print("--- Results ---")
    print("Expected revenue for each protection level (k):")
    for k_val, rev in all_k_revenues.items():
        print(f"k = {k_val:2d}: ${rev:8.2f}")

    print(f"\nMaximum expected revenue is ${max_revenue:.2f} with a protection level of k = {best_k}.\n")

    # --- Part 2: Calculate ticket sale details for the optimal k ---
    print(f"Step 2: Analyzing the optimal policy (k={best_k}) to form the final equation...")
    
    # DP table storing (revenue, expected_cheap_tickets, expected_expensive_tickets)
    Ed = [[(0.0, 0.0, 0.0)] * (SEAT_CAPACITY + 1) for _ in range(NUM_DAYS + 1)]

    for t in range(1, NUM_DAYS + 1):
        day_index = NUM_DAYS - t + 1
        is_second_week = day_index > 7

        @functools.lru_cache(maxsize=None)
        def get_daily_detailed_stats(n1, n2, s_avail):
            if s_avail <= 0 or (n1 == 0 and n2 == 0):
                return Ed[t - 1][max(0, s_avail)]

            if n2 > 0:
                if s_avail > best_k:
                    rev, c, e = get_daily_detailed_stats(n1, n2 - 1, s_avail - 1)
                    return PRICE_CHEAP + rev, 1 + c, e
                else:
                    rev_b, c_b, e_b = get_daily_detailed_stats(n1, n2 - 1, s_avail - 1)
                    rev_w, c_w, e_w = get_daily_detailed_stats(n1, n2 - 1, s_avail)
                    
                    exp_rev = PROB_BUY_EXPENSIVE * (PRICE_EXPENSIVE + rev_b) + (1 - PROB_BUY_EXPENSIVE) * rev_w
                    exp_c = PROB_BUY_EXPENSIVE * c_b + (1 - PROB_BUY_EXPENSIVE) * c_w
                    exp_e = PROB_BUY_EXPENSIVE * (1 + e_b) + (1 - PROB_BUY_EXPENSIVE) * e_w
                    return exp_rev, exp_c, exp_e
            
            if n1 > 0:
                if s_avail > best_k:
                    rev, c, e = get_daily_detailed_stats(n1 - 1, n2, s_avail - 1)
                    return PRICE_CHEAP + rev, 1 + c, e
                else:
                    return get_daily_detailed_stats(n1 - 1, n2, s_avail)
            return Ed[t-1][s_avail] # Should not be reached if n1/n2>0

        for s in range(1, SEAT_CAPACITY + 1):
            exp_rev, exp_c, exp_e = 0.0, 0.0, 0.0
            p2_dist_today = P2_DIST if is_second_week else {0: 1.0}

            for n1_arrivals, p1 in P1_DIST.items():
                for n2_arrivals, p2 in p2_dist_today.items():
                    prob = p1 * p2
                    r, c, e = get_daily_detailed_stats(n1_arrivals, n2_arrivals, s)
                    exp_rev += prob * r
                    exp_c += prob * c
                    exp_e += prob * e
            
            Ed[t][s] = (exp_rev, exp_c, exp_e)
        
        get_daily_detailed_stats.cache_clear()

    final_stats = Ed[NUM_DAYS][SEAT_CAPACITY]
    expected_cheap_sold = final_stats[1]
    expected_expensive_sold = final_stats[2]
    
    print("\nThe maximum expected revenue can be broken down as follows:")
    print("Final Equation based on expected ticket sales:")
    print(f"  Expected cheap tickets sold: {expected_cheap_sold:.4f}")
    print(f"  Expected expensive tickets sold: {expected_expensive_sold:.4f}")
    print(f"  {expected_cheap_sold:.4f} * {PRICE_CHEAP} + {expected_expensive_sold:.4f} * {PRICE_EXPENSIVE} = ${max_revenue:.2f}")

    return max_revenue

if __name__ == '__main__':
    max_rev = solve_flight_revenue()
    print(f"\nFinal Answer: The maximum expected total revenue is ${max_rev:.2f}.")
    print(f"<<<{max_rev:.2f}>>>")
