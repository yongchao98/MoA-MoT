import math

def solve_revenue_maximization():
    """
    This function calculates the maximum expected revenue for an airline
    using dynamic programming.

    The state of the DP is (t, s), where:
    - t is the number of days remaining until departure (from 1 to 14).
    - s is the number of seats available (from 0 to 10).

    V[t][s] stores the maximum expected revenue from day t onwards with s seats.
    policy[t][s] stores the optimal decision (number of cheap tickets to offer).

    The calculation proceeds backwards from t=1 to t=14.
    """

    # DP table: V[t][s]
    V = [[0.0 for _ in range(11)] for _ in range(15)]
    # Optimal policy table: policy[t][s]
    policy = [[0 for _ in range(11)] for _ in range(15)]

    # Demand distributions
    P_d1 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1 demand
    P_d2 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2 demand (week 2)
    P_d2_none = {0: 1.0}               # Class 2 demand (week 1)

    # Helper for combinations
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    # Main DP loop
    for t in range(1, 15):  # t = days left, from 1 to 14
        # Class 2 customers only appear in the last 7 days (t <= 7)
        d2_dist = P_d2 if t <= 7 else P_d2_none

        for s in range(1, 11):  # s = seats remaining
            max_rev_for_s = -1.0
            best_c = -1

            # Iterate through all possible decisions `c` (number of cheap tickets to offer)
            for c in range(s + 1):
                e_avail = s - c  # Number of expensive tickets available
                current_total_expected_rev = 0.0

                # Iterate through all possible demand scenarios for d1 and d2
                for d1, prob_d1 in P_d1.items():
                    for d2, prob_d2 in d2_dist.items():
                        
                        prob_scenario = prob_d1 * prob_d2
                        
                        # Revenue from cheap tickets sold today
                        sold_c_today = min(c, d1 + d2)
                        rev_c_today = 100 * sold_c_today
                        
                        # Class 2 customers who couldn't get a cheap ticket
                        d2_rem = max(0, d2 - c)
                        
                        # Calculate expected revenue from expensive sales and future value
                        exp_rev_from_exp_and_future = 0.0
                        
                        if d2_rem > 0 and e_avail > 0:
                            # Case: Potential for expensive ticket sales
                            # Iterate over the number of C2 customers (j) who decide to buy an expensive ticket (prob 0.5)
                            for j in range(d2_rem + 1):
                                prob_j = combinations(d2_rem, j) * (0.5 ** d2_rem)
                                
                                sold_e_today = min(j, e_avail)
                                rev_e_today = 200 * sold_e_today
                                
                                seats_left_end_of_day = s - sold_c_today - sold_e_today
                                future_rev = V[t-1][seats_left_end_of_day]
                                
                                exp_rev_from_exp_and_future += prob_j * (rev_e_today + future_rev)
                        else:
                            # Case: No potential for expensive sales today
                            seats_left_end_of_day = s - sold_c_today
                            future_rev = V[t-1][seats_left_end_of_day]
                            exp_rev_from_exp_and_future = future_rev

                        current_total_expected_rev += prob_scenario * (rev_c_today + exp_rev_from_exp_and_future)

                if current_total_expected_rev > max_rev_for_s:
                    max_rev_for_s = current_total_expected_rev
                    best_c = c
            
            V[t][s] = max_rev_for_s
            policy[t][s] = best_c

    # --- Outputting the result ---
    
    # The final answer is V[14][10]
    final_answer = V[14][10]
    
    # Find the optimal decision on the first day (t=14, s=10)
    c_star = policy[14][10]

    # Calculate the components of the expectation for the final equation
    # On day 14, d2=0, so we only consider d1.
    # term_k = (revenue from d1=k) + (expected future revenue)
    term0 = 100 * min(c_star, 0) + V[13][10 - min(c_star, 0)]
    term1 = 100 * min(c_star, 1) + V[13][10 - min(c_star, 1)]
    term2 = 100 * min(c_star, 2) + V[13][10 - min(c_star, 2)]

    print("The maximum expected total revenue is found using dynamic programming.")
    print("The optimal decision on the first day (t=14, 10 seats) is to offer a certain number of cheap tickets.")
    print(f"The optimal number of cheap tickets to offer is {c_star}.")
    print("\nThe expected revenue for this optimal decision is calculated by considering the demand distribution for Class 1 customers (P(0)=0.25, P(1)=0.5, P(2)=0.25):")
    print("E[Revenue] = 0.25 * E[Rev|0 requests] + 0.50 * E[Rev|1 request] + 0.25 * E[Rev|2 requests]")
    print("\nThe final equation with the calculated values is:")
    print(f"0.25 * {term0:.2f} + 0.50 * {term1:.2f} + 0.25 * {term2:.2f} = {final_answer:.2f}")
    
    print(f"<<<{final_answer:.2f}>>>")

solve_revenue_maximization()