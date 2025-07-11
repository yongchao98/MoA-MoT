import math

def solve_revenue_maximization():
    """
    This function calculates the maximum expected revenue for an airline by optimizing its booking policy.
    The policy is static: we set a fixed number of cheap tickets (q) and expensive tickets (10-q)
    to be available over the 14-day sales period. We iterate through all possible values of q (0 to 10)
    to find the one that maximizes the expected revenue.

    For each policy q, we use dynamic programming to calculate the total expected revenue.
    The state of our DP is defined by V[d][c][e], representing the expected revenue from day 'd' onwards,
    given that 'c' cheap and 'e' expensive tickets are still available.

    The DP works backward from the last day (14) to the first (1). On each day, we consider all
    possible customer arrival scenarios and their probabilities to calculate the expected revenue
    for that day plus the expected future revenue from the next state.
    """

    max_total_revenue = -1.0
    best_q = -1
    best_exp_cheap_sales = 0.0
    best_exp_exp_sales = 0.0

    p1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    p2_dist_week2 = {0: 0.25, 1: 0.5, 2: 0.25}
    p2_dist_week1 = {0: 1.0}  # Only 0 arrivals with prob 1.0

    # Iterate through all possible static policies (q = number of cheap tickets)
    for q in range(11):
        e_cap = 10 - q
        
        # DP tables for this policy q
        # V[d][c][e]: expected revenue from day d onwards
        # E_c[d][c][e]: expected number of cheap tickets sold from day d onwards
        # E_e[d][c][e]: expected number of expensive tickets sold from day d onwards
        V = [[[0.0 for _ in range(e_cap + 1)] for _ in range(q + 1)] for _ in range(16)]
        E_c = [[[0.0 for _ in range(e_cap + 1)] for _ in range(q + 1)] for _ in range(16)]
        E_e = [[[0.0 for _ in range(e_cap + 1)] for _ in range(q + 1)] for _ in range(16)]
        
        # Iterate days backwards from d=14 to 1
        for d in range(14, 0, -1):
            p2_dist = p2_dist_week2 if d >= 8 else p2_dist_week1
            
            for c in range(q + 1):
                for e in range(e_cap + 1):
                    
                    exp_val_for_state = 0.0
                    exp_c_for_state = 0.0
                    exp_e_for_state = 0.0
                    
                    # Iterate over demand scenarios for Class 1 (r1) and Class 2 (r2)
                    for r1, p1 in p1_dist.items():
                        for r2, p2 in p2_dist.items():
                            prob_scenario = p1 * p2
                            
                            # Class 2 has priority for cheap tickets
                            sold_c2 = min(r2, c)
                            c_after_c2 = c - sold_c2
                            r2_rem = r2 - sold_c2
                            
                            sold_c1 = min(r1, c_after_c2)
                            c_final = c_after_c2 - sold_c1
                            
                            daily_c_sold = sold_c1 + sold_c2
                            daily_c_rev = daily_c_sold * 100

                            # Handle probabilistic expensive sales
                            exp_val_from_exp_branch = 0.0
                            exp_c_from_exp_branch = 0.0
                            exp_e_from_exp_branch = 0.0

                            num_c2_consider_exp = r2_rem
                            
                            # Iterate over k, the number of C2 customers who decide to buy expensive
                            for k in range(num_c2_consider_exp + 1):
                                prob_k = math.comb(num_c2_consider_exp, k) * (0.5 ** num_c2_consider_exp)
                                
                                sold_e2 = min(k, e)
                                e_final = e - sold_e2
                                daily_e_rev = sold_e2 * 200

                                path_val = daily_c_rev + daily_e_rev + V[d+1][c_final][e_final]
                                exp_val_from_exp_branch += prob_k * path_val
                                
                                path_c_sold = daily_c_sold + E_c[d+1][c_final][e_final]
                                exp_c_from_exp_branch += prob_k * path_c_sold
                                
                                path_e_sold = sold_e2 + E_e[d+1][c_final][e_final]
                                exp_e_from_exp_branch += prob_k * path_e_sold

                            exp_val_for_state += prob_scenario * exp_val_from_exp_branch
                            exp_c_for_state += prob_scenario * exp_c_from_exp_branch
                            exp_e_for_state += prob_scenario * exp_e_from_exp_branch

                    V[d][c][e] = exp_val_for_state
                    E_c[d][c][e] = exp_c_for_state
                    E_e[d][c][e] = exp_e_for_state
        
        total_revenue_for_q = V[1][q][e_cap]
        if total_revenue_for_q > max_total_revenue:
            max_total_revenue = total_revenue_for_q
            best_q = q
            best_exp_cheap_sales = E_c[1][q][e_cap]
            best_exp_exp_sales = E_e[1][q][e_cap]
    
    # Print the final result and the equation
    print(f"The maximum expected revenue is {max_total_revenue:.2f}.")
    print(f"This is achieved by setting the total number of cheap tickets to {best_q}.")
    print("\nThe final revenue calculation is:")
    print(f"{best_exp_cheap_sales:.4f} * 100 + {best_exp_exp_sales:.4f} * 200 = {max_total_revenue:.2f}")

solve_revenue_maximization()