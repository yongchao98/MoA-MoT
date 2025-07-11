import math

def solve_flight_revenue():
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    
    memo_comb = {}
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        if (n, k) in memo_comb:
            return memo_comb[(n, k)]
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        memo_comb[(n, k)] = res
        return res

    # DP table: E[t][c][e]
    # t: day (1 to 15, 1-indexed)
    # c: remaining cheap tickets (0 to 10)
    # e: remaining expensive tickets (0 to 10)
    E = [[[0.0 for _ in range(11)] for _ in range(11)] for _ in range(16)]

    # Demand distributions
    d1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    d2_dist_week1 = {0: 1.0}
    d2_dist_week2 = {0: 0.25, 1: 0.5, 2: 0.25}

    # DP calculation, iterating backwards in time from day 14 to 1
    for t in range(14, 0, -1):
        current_d2_dist = d2_dist_week1 if t <= 7 else d2_dist_week2

        for c in range(11):  # Remaining cheap tickets
            for e in range(11):  # Remaining expensive tickets
                
                expected_value_for_state = 0.0
                
                # Iterate over all possible demand scenarios for the day
                for d1, p1 in d1_dist.items():
                    for d2, p2 in current_d2_dist.items():
                        prob_scenario = p1 * p2
                        
                        # --- Simulate sales for the (d1, d2) scenario ---
                        
                        # 1. Class 2 customers (priority) purchase cheap tickets
                        s2_cheap = min(d2, c)
                        rev_c2_cheap = s2_cheap * 100
                        c_after_c2 = c - s2_cheap
                        d2_rem_for_exp = d2 - s2_cheap
                        
                        # 2. Class 1 customers purchase cheap tickets
                        s1_cheap = min(d1, c_after_c2)
                        rev_c1_cheap = s1_cheap * 100
                        c_final = c_after_c2 - s1_cheap
                        
                        immediate_cheap_rev = rev_c1_cheap + rev_c2_cheap
                        
                        # 3. Class 2 customers consider expensive tickets
                        n_potential_buyers = min(d2_rem_for_exp, e)
                        exp_val_from_exp_sales_and_future = 0.0
                        
                        # Calculate expected value over the stochastic expensive sales
                        for k in range(n_potential_buyers + 1):
                            prob_k_sales = combinations(n_potential_buyers, k) * (0.5)**n_potential_buyers
                            
                            immediate_exp_rev = k * 200
                            e_final = e - k
                            
                            future_rev = E[t + 1][c_final][e_final]
                            
                            exp_val_from_exp_sales_and_future += prob_k_sales * (immediate_exp_rev + future_rev)

                        scenario_value = immediate_cheap_rev + exp_val_from_exp_sales_and_future
                        expected_value_for_state += prob_scenario * scenario_value

                E[t][c][e] = expected_value_for_state

    # Find the optimal initial allocation of 10 seats
    max_revenue = 0.0
    optimal_policy = (0, 0)
    total_seats = 10
    
    for c_init in range(total_seats + 1):
        e_init = total_seats - c_init
        current_revenue = E[1][c_init][e_init]
        
        if current_revenue > max_revenue:
            max_revenue = current_revenue
            optimal_policy = (c_init, e_init)
            
    opt_c, opt_e = optimal_policy
    
    print(f"The optimal policy is to offer {opt_c} cheap tickets and {opt_e} expensive tickets.")
    print(f"The maximum expected total revenue is: {max_revenue:.2f}")

solve_flight_revenue()