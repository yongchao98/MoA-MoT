import math

def solve_flight_revenue():
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    CAPACITY = 10
    DAYS = 14
    CHEAP_PRICE = 100
    EXPENSIVE_PRICE = 200

    # Customer arrival probabilities: P(0)=0.25, P(1)=0.5, P(2)=0.25
    PROBS = {0: 0.25, 1: 0.5, 2: 0.25}

    # dp[d][s]: max expected revenue with d days left and s seats available.
    dp = [[0.0] * (CAPACITY + 1) for _ in range(DAYS + 1)]
    # policy[d][s]: optimal number of cheap tickets to offer for state (d,s).
    policy = [[0] * (CAPACITY + 1) for _ in range(DAYS + 1)]

    # --- Dynamic Programming Calculation ---
    # Iterate from Day 1 (last day of sales) to Day 14 (first day).
    for d in range(1, DAYS + 1):
        for s in range(CAPACITY + 1):
            max_rev_for_state = -1.0
            best_c = -1

            # Iterate through all possible policies: 'c' is the number of cheap tickets to offer.
            for c in range(s + 1):
                expected_rev_for_policy_c = 0.0

                # Determine customer arrival distributions for the current day.
                # Class 2 customers only arrive in the second week (d <= 7).
                is_week_one = d > 7
                n1_arrivals = PROBS.keys()
                n2_arrivals = [0] if is_week_one else PROBS.keys()

                # Sum expected revenue over all arrival scenarios (n1, n2).
                for n1 in n1_arrivals:
                    for n2 in n2_arrivals:
                        prob_scenario = PROBS[n1] * (1 if is_week_one else PROBS[n2])

                        # 1. Cheap ticket sales (Class 2 has priority).
                        s2_cheap = min(n2, c)
                        s1_cheap = min(n1, c - s2_cheap)
                        cheap_revenue = CHEAP_PRICE * (s1_cheap + s2_cheap)
                        
                        seats_sold_cheap = s1_cheap + s2_cheap
                        seats_after_cheap = s - seats_sold_cheap
                        n2_rem_for_exp = n2 - s2_cheap

                        # 2. Expensive ticket sales (only Class 2).
                        # Calculate the expected value from this phase, considering the 50% purchase probability.
                        expected_value_from_exp_phase = 0.0
                        if n2_rem_for_exp > 0 and seats_after_cheap > 0:
                            # 'k' is the number of customers (out of n2_rem_for_exp) who decide to buy expensive.
                            for k in range(n2_rem_for_exp + 1):
                                # Binomial probability: P(k successes in n trials with p=0.5)
                                prob_k = math.comb(n2_rem_for_exp, k) * (0.5 ** n2_rem_for_exp)
                                
                                s2_exp = min(k, seats_after_cheap)
                                expensive_revenue = EXPENSIVE_PRICE * s2_exp
                                seats_final = seats_after_cheap - s2_exp
                                future_revenue = dp[d - 1][seats_final]
                                
                                expected_value_from_exp_phase += prob_k * (expensive_revenue + future_revenue)
                        else:
                            # No expensive sales if no remaining Class 2 customers or no seats.
                            seats_final = seats_after_cheap
                            future_revenue = dp[d - 1][seats_final]
                            expected_value_from_exp_phase = future_revenue

                        # Total expected value for this (n1, n2) scenario.
                        total_value_for_scenario = cheap_revenue + expected_value_from_exp_phase
                        expected_rev_for_policy_c += prob_scenario * total_value_for_scenario

                # Check if this policy 'c' is the best so far for state (d,s).
                if expected_rev_for_policy_c > max_rev_for_state:
                    max_rev_for_state = expected_rev_for_policy_c
                    best_c = c
            
            dp[d][s] = max_rev_for_state
            policy[d][s] = best_c

    # --- Final Result Presentation ---
    # The final answer is the max expected revenue at the start: Day 14 with 10 seats.
    total_revenue = dp[DAYS][CAPACITY]
    
    # Deconstruct the total revenue for the first day (d=14, s=10) to show the "equation".
    d_start = DAYS
    s_start = CAPACITY
    c_optimal_start = policy[d_start][s_start]
    
    exp_rev_day_14 = 0.0
    exp_future_rev = 0.0

    # On Day 14 (first week), only Class 1 customers arrive.
    for n1 in PROBS.keys():
        prob_n1 = PROBS[n1]
        
        # Sales on Day 14 are simple: only Class 1 buying cheap tickets.
        s1_cheap = min(n1, c_optimal_start)
        day_revenue = CHEAP_PRICE * s1_cheap
        
        seats_final = s_start - s1_cheap
        future_revenue = dp[d_start - 1][seats_final]

        exp_rev_day_14 += prob_n1 * day_revenue
        exp_future_rev += prob_n1 * future_revenue
    
    # Print the final result and the breakdown.
    print(f"The Maximum Expected Total Revenue is {total_revenue:.2f}.")
    print(f"The optimal policy for the first day of sales (Day 14) with {CAPACITY} seats is to offer {c_optimal_start} cheap tickets.")
    print("\nThe final revenue calculation breaks down as follows:")
    # Using the numbers calculated from the optimal policy on the first day
    first_day_rev = exp_rev_day_14
    all_other_days_rev = exp_future_rev
    print(f"{first_day_rev:.2f} + {all_other_days_rev:.2f} = {first_day_rev + all_other_days_rev:.2f}")

solve_flight_revenue()
<<<1114.35>>>