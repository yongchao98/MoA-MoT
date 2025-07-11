import collections

# Memoization for the sub-problem of calculating expensive purchase distribution
memo_exp_purchase = {}

def get_exp_purchase_dist(num_customers, num_seats):
    """
    Calculates the probability distribution of the number of expensive tickets sold
    to a group of customers who decide with 50% probability.

    Args:
        num_customers: Number of Class 2 customers who were denied a cheap ticket.
        num_seats: Number of available seats for them to purchase as expensive.

    Returns:
        A dictionary {num_sales: probability}.
    """
    # Use memoization to avoid re-computing results
    if (num_customers, num_seats) in memo_exp_purchase:
        return memo_exp_purchase[(num_customers, num_seats)]

    # Base cases: no customers or no seats means no sales.
    if num_customers == 0 or num_seats == 0:
        return {0: 1.0}

    # Recursive step: consider the first customer's decision.
    # Case 1: Customer walks away (probability 0.5).
    # The remaining (num_customers - 1) customers face the same `num_seats`.
    dist_walk_away = get_exp_purchase_dist(num_customers - 1, num_seats)

    # Case 2: Customer tries to buy (probability 0.5) and succeeds.
    # The remaining (num_customers - 1) customers face (num_seats - 1).
    dist_buy_raw = get_exp_purchase_dist(num_customers - 1, num_seats - 1)
    # Add 1 to each sale count in the resulting distribution.
    dist_buy = {sales + 1: prob for sales, prob in dist_buy_raw.items()}

    # Combine the two distributions, weighted by their 0.5 probability.
    final_dist = collections.defaultdict(float)
    for sales, prob in dist_walk_away.items():
        final_dist[sales] += 0.5 * prob
    for sales, prob in dist_buy.items():
        final_dist[sales] += 0.5 * prob
        
    memo_exp_purchase[(num_customers, num_seats)] = dict(final_dist)
    return memo_exp_purchase[(num_customers, num_seats)]

def solve():
    """
    Main function to solve the revenue management problem.
    """
    max_revenue = -1.0
    best_k = -1
    
    # Customer demand distributions
    D1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    D2_dist_week2 = {0: 0.25, 1: 0.5, 2: 0.25}
    D2_dist_week1 = {0: 1.0}

    # Iterate through all possible protection levels (k)
    # k = number of seats protected for expensive fares
    for k in range(11):  # k from 0 to 10
        cheap_quota = 10 - k
        
        # DP table: V[d, s_c, s_e] = max expected revenue from day d onwards
        V = collections.defaultdict(float)

        # Loop backwards in time, from day 1 to day 14
        for d in range(1, 15):
            current_D2_dist = D2_dist_week2 if d <= 7 else D2_dist_week1

            for s_c in range(cheap_quota + 1):
                for s_e in range(11):
                    if s_c > s_e: continue

                    expected_value_for_state = 0
                    for d1, p1 in D1_dist.items():
                        for d2, p2 in current_D2_dist.items():
                            prob_scenario = p1 * p2
                            
                            s_c_initial, s_e_initial = s_c, s_e
                            c2_sold_cheap = min(d2, s_c_initial, s_e_initial)
                            rev_c2_cheap = c2_sold_cheap * 100
                            s_c_after_c2, s_e_after_c2 = s_c_initial - c2_sold_cheap, s_e_initial - c2_sold_cheap
                            d2_rem = d2 - c2_sold_cheap

                            c1_sold_cheap = min(d1, s_c_after_c2, s_e_after_c2)
                            rev_c1_cheap = c1_sold_cheap * 100
                            s_c_after_c1, s_e_after_c1 = s_c_after_c2 - c1_sold_cheap, s_e_after_c2 - c1_sold_cheap
                            
                            exp_purchase_dist = get_exp_purchase_dist(d2_rem, s_e_after_c1)
                            
                            path_value = 0
                            for c2_sold_exp, prob_exp in exp_purchase_dist.items():
                                rev_c2_exp = c2_sold_exp * 200
                                s_c_final, s_e_final = s_c_after_c1, s_e_after_c1 - c2_sold_exp
                                future_rev = V[d - 1, s_c_final, s_e_final]
                                total_rev_for_outcome = rev_c1_cheap + rev_c2_cheap + rev_c2_exp + future_rev
                                path_value += prob_exp * total_rev_for_outcome
                            
                            expected_value_for_state += prob_scenario * path_value
                    
                    V[d, s_c, s_e] = expected_value_for_state

        current_revenue = V[14, cheap_quota, 10]
        if current_revenue > max_revenue:
            max_revenue = current_revenue
            best_k = k

    # Forward pass to calculate revenue components for the optimal policy
    cheap_quota = 10 - best_k
    State_Prob = collections.defaultdict(float)
    State_Prob[14, cheap_quota, 10] = 1.0
    
    E_rev_c1_cheap, E_rev_c2_cheap, E_rev_c2_exp = 0, 0, 0

    for d in range(14, 0, -1):
        current_D2_dist = D2_dist_week2 if d <= 7 else D2_dist_week1
        Next_State_Prob = collections.defaultdict(float)

        for (s_c, s_e), prob_state in State_Prob.items():
            if prob_state == 0: continue

            for d1, p1 in D1_dist.items():
                for d2, p2 in current_D2_dist.items():
                    prob_scenario = p1 * p2
                    
                    s_c_initial, s_e_initial = s_c, s_e
                    c2_sold_cheap = min(d2, s_c_initial, s_e_initial)
                    s_c_after_c2, s_e_after_c2 = s_c_initial - c2_sold_cheap, s_e_initial - c2_sold_cheap
                    d2_rem = d2 - c2_sold_cheap

                    c1_sold_cheap = min(d1, s_c_after_c2, s_e_after_c2)
                    s_c_after_c1, s_e_after_c1 = s_c_after_c2 - c1_sold_cheap, s_e_after_c2 - c1_sold_cheap
                    
                    E_rev_c1_cheap += prob_state * prob_scenario * (c1_sold_cheap * 100)
                    E_rev_c2_cheap += prob_state * prob_scenario * (c2_sold_cheap * 100)
                    
                    exp_purchase_dist = get_exp_purchase_dist(d2_rem, s_e_after_c1)
                    
                    for c2_sold_exp, prob_exp in exp_purchase_dist.items():
                        E_rev_c2_exp += prob_state * prob_scenario * prob_exp * (c2_sold_exp * 200)
                        s_c_final, s_e_final = s_c_after_c1, s_e_after_c1 - c2_sold_exp
                        Next_State_Prob[s_c_final, s_e_final] += prob_state * prob_scenario * prob_exp
        
        State_Prob = Next_State_Prob

    print(f"The maximum expected total revenue is achieved by protecting {best_k} seats for expensive fares (i.e., setting a quota of {10 - best_k} cheap tickets).")
    print("The breakdown of the maximum expected revenue is as follows:")
    print(f"Expected Revenue from Class 1 (cheap tickets) = {E_rev_c1_cheap:.2f}")
    print(f"Expected Revenue from Class 2 (cheap tickets) = {E_rev_c2_cheap:.2f}")
    print(f"Expected Revenue from Class 2 (expensive tickets) = {E_rev_c2_exp:.2f}")
    print("The final equation for the maximum expected revenue is:")
    print(f"{E_rev_c1_cheap:.2f} + {E_rev_c2_cheap:.2f} + {E_rev_c2_exp:.2f} = {max_revenue:.2f}")
    
    print(f"<<<{max_revenue:.2f}>>>")

solve()