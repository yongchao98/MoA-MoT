import math

def solve_revenue_management():
    """
    Calculates the maximum expected total revenue using dynamic programming.
    """
    NUM_DAYS = 14
    CAPACITY = 10
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200

    # Demand probabilities
    # P(D=k) = P_dist[k]
    P1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    P2_active_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    P2_inactive_dist = {0: 1.0} # No class 2 customers

    # DP table: V[t][s] = max expected revenue with t days left and s seats
    # We use t as days remaining, so t=1 is the last day of sale.
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(NUM_DAYS + 1)]

    # Base case V[0][s] = 0 is already initialized.

    # Iterate backwards in time (from t=1 day remaining up to t=14 days remaining)
    for t in range(1, NUM_DAYS + 1):
        # Determine Class 2 customer demand for this day
        # Last 7 days (t=1 to 7) have active Class 2 customers
        if t <= 7:
            P2_dist = P2_active_dist
        else: # First 7 days (t=8 to 14) have no Class 2 customers
            P2_dist = P2_inactive_dist
        
        # Iterate over number of available seats
        for s in range(1, CAPACITY + 1):
            max_exp_rev_for_s = 0.0
            
            # Iterate over the decision variable p (protection level)
            for p in range(s + 1):
                cheap_tickets_available = s - p
                current_exp_rev_for_p = 0.0

                # Iterate over all possible demand scenarios for Class 1 (d1) and Class 2 (d2)
                for d1, prob1 in P1_dist.items():
                    for d2, prob2 in P2_dist.items():
                        prob_d1d2 = prob1 * prob2
                        
                        # --- Simulate sales for this (d1, d2) scenario ---
                        
                        # 1. Class 2 customers buy cheap tickets (priority)
                        sold_c2 = min(d2, cheap_tickets_available)
                        rev_cheap = PRICE_CHEAP * sold_c2
                        
                        s_after_c2 = s - sold_c2
                        c_after_c2 = cheap_tickets_available - sold_c2
                        d2_unfulfilled = d2 - sold_c2
                        
                        # 2. Class 1 customers buy remaining cheap tickets
                        sold_c1 = min(d1, c_after_c2)
                        rev_cheap += PRICE_CHEAP * sold_c1
                        
                        s_after_c1 = s_after_c2 - sold_c1
                        
                        # 3. Class 2 customers (who couldn't get cheap) consider expensive tickets
                        n_exp_candidates = min(d2_unfulfilled, s_after_c1)
                        
                        # This part calculates E[rev_expensive + V(t-1, s_next)]
                        # over the binomial outcomes of expensive sales
                        exp_value_from_exp_sales = 0.0
                        
                        if n_exp_candidates > 0:
                            prob_k_factor = 0.5 ** n_exp_candidates
                            for k in range(n_exp_candidates + 1):
                                # k is the number of expensive tickets sold
                                prob_k = math.comb(n_exp_candidates, k) * prob_k_factor
                                rev_expensive = PRICE_EXPENSIVE * k
                                s_next = s_after_c1 - k
                                future_rev = V[t - 1][s_next]
                                exp_value_from_exp_sales += prob_k * (rev_expensive + future_rev)
                        else:
                            # If no candidates, no expensive sales. Future revenue is based on current state.
                            s_next = s_after_c1
                            future_rev = V[t - 1][s_next]
                            exp_value_from_exp_sales = future_rev

                        # Total expected value for this (d1, d2) pair
                        current_exp_rev_for_p += prob_d1d2 * (rev_cheap + exp_value_from_exp_sales)

                # Update max revenue for state (t,s)
                if current_exp_rev_for_p > max_exp_rev_for_s:
                    max_exp_rev_for_s = current_exp_rev_for_p
            
            V[t][s] = max_exp_rev_for_s

    # The final answer is the max expected revenue starting with 14 days and 10 seats
    final_answer = V[NUM_DAYS][CAPACITY]
    print(f"The maximum expected total revenue is: {final_answer:.2f}")

if __name__ == '__main__':
    solve_revenue_management()
