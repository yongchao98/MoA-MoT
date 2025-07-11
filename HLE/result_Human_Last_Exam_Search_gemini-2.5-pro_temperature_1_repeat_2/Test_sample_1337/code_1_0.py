import math

def nCr(n, r):
    """Calculates the binomial coefficient 'n choose r'."""
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def solve_revenue_management():
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    # V[d][s]: max expected revenue from day d to 14, with s seats.
    # d is 1-indexed (1 to 15). Day 15 is after the selling period.
    # s is the number of seats, 0 to 10.
    V = [[0.0] * 11 for _ in range(16)]
    
    # P1_dist: Probability distribution for Class 1 customer arrivals.
    P1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    # P2_dist: Probability distribution for Class 2 customer arrivals.
    P2_dist = {0: 0.25, 1: 0.5, 2: 0.25}

    # Iterate backwards from the last day of sales (d=14) to the first (d=1).
    for d in range(14, 0, -1):
        # Iterate over the number of available seats.
        for s in range(0, 11):
            max_rev_for_s = -1.0

            # Iterate over all possible protection levels p.
            for p in range(0, s + 1):
                current_p_exp_rev = 0.0
                
                # Determine Class 2 customer distribution for the current day.
                # Class 2 customers only appear in the second week (days 8-14).
                is_week_2 = d > 7
                d2_distribution = P2_dist if is_week_2 else {0: 1.0}

                # Iterate over all possible demand scenarios for Class 1.
                for d1, prob1 in P1_dist.items():
                    # Iterate over all possible demand scenarios for Class 2.
                    for d2, prob2 in d2_distribution.items():
                        
                        prob_scenario = prob1 * prob2
                        
                        cheap_tickets_available = s - p
                        
                        # --- Cheap Ticket Sales ---
                        sold_cheap_c2 = min(d2, cheap_tickets_available)
                        remaining_cheap_for_c1 = cheap_tickets_available - sold_cheap_c2
                        sold_cheap_c1 = min(d1, remaining_cheap_for_c1)
                        total_sold_cheap = sold_cheap_c1 + sold_cheap_c2
                        revenue_from_cheap = 100 * total_sold_cheap

                        # --- Expensive Ticket Sales ---
                        unaccommodated_c2 = d2 - sold_cheap_c2
                        potential_expensive_buyers = min(unaccommodated_c2, p)
                        
                        exp_value_from_expensive_sales_and_future = 0.0
                        
                        if potential_expensive_buyers > 0:
                            prob_binom_base = 0.5 ** potential_expensive_buyers
                            for k in range(potential_expensive_buyers + 1):
                                prob_k = nCr(potential_expensive_buyers, k) * prob_binom_base
                                revenue_from_expensive = 200 * k
                                seats_remaining_for_next_day = s - total_sold_cheap - k
                                
                                future_revenue = V[d + 1][seats_remaining_for_next_day] if seats_remaining_for_next_day >= 0 else 0.0
                                exp_value_from_expensive_sales_and_future += prob_k * (revenue_from_expensive + future_revenue)
                        else:
                            seats_remaining_for_next_day = s - total_sold_cheap
                            future_revenue = V[d + 1][seats_remaining_for_next_day] if seats_remaining_for_next_day >= 0 else 0.0
                            exp_value_from_expensive_sales_and_future = future_revenue
                            
                        current_p_exp_rev += prob_scenario * (revenue_from_cheap + exp_value_from_expensive_sales_and_future)
                
                if current_p_exp_rev > max_rev_for_s:
                    max_rev_for_s = current_p_exp_rev
            
            V[d][s] = max_rev_for_s
            
    return V[1][10]

# Calculate the final answer.
# max_revenue = solve_revenue_management()
# print(f"{max_revenue:.2f}")