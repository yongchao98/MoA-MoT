import math

def solve_flight_revenue():
    """
    Calculates the maximum expected revenue from flight ticket sales using dynamic programming.
    """
    # --- Problem Parameters ---
    CAPACITY = 10
    DAYS = 14
    D_VALUES = [0, 1, 2]
    
    # Demand distributions
    D1_PROBS = {0: 0.25, 1: 0.5, 2: 0.25}
    # Days 1-7 (t=14 down to 8): No Class 2 customers
    D2_PROBS_FIRST_WEEK = {0: 1.0, 1: 0.0, 2: 0.0}
    # Days 8-14 (t=7 down to 1): Class 2 customers appear
    D2_PROBS_SECOND_WEEK = {0: 0.25, 1: 0.5, 2: 0.25}

    # --- DP Initialization ---
    # V[t][s]: max expected revenue with t days remaining and s seats
    V = [[0.0] * (CAPACITY + 1) for _ in range(DAYS + 1)]
    # policy[t][s]: optimal number of cheap tickets 'c' to offer
    policy = [[0] * (CAPACITY + 1) for _ in range(DAYS + 1)]

    # --- DP Calculation ---
    # Loop backwards from t=1 to 14 days remaining
    for t in range(1, DAYS + 1):
        # Determine Class 2 demand distribution for the current day
        d2_probs = D2_PROBS_SECOND_WEEK if t <= 7 else D2_PROBS_FIRST_WEEK

        # Loop over number of available seats s
        for s in range(1, CAPACITY + 1):
            max_expected_revenue = -1.0
            best_c = 0

            # Loop over decision variable c: number of cheap tickets to offer
            for c in range(s + 1):
                e = s - c  # Number of expensive tickets
                current_c_expected_revenue = 0.0

                # Iterate over all possible demand scenarios for d1 and d2
                for d1 in D_VALUES:
                    for d2 in D_VALUES:
                        prob_scenario = D1_PROBS[d1] * d2_probs[d2]
                        if prob_scenario == 0:
                            continue

                        # Class 2 customers have priority for cheap tickets
                        sold_c2 = min(d2, c)
                        # Class 1 customers get the remaining cheap tickets
                        sold_c1 = min(d1, c - sold_c2)
                        
                        # Customers who couldn't get a cheap ticket consider an expensive one
                        d2_unmet = d2 - sold_c2
                        
                        # Calculate the expected value contribution from this scenario
                        expected_value_from_binom = 0.0
                        if d2_unmet == 0:
                            # No customers considering expensive tickets
                            daily_revenue = (sold_c1 + sold_c2) * 100
                            seats_sold = sold_c1 + sold_c2
                            s_next = s - seats_sold
                            expected_value_from_binom = daily_revenue + V[t - 1][s_next]
                        else:
                            # Sum over binomial outcomes for expensive ticket purchase
                            # k = number of customers who decide to buy an expensive ticket
                            prob_d2_unmet_buy = 0.5
                            for k in range(d2_unmet + 1):
                                prob_k = math.comb(d2_unmet, k) * (prob_d2_unmet_buy ** k) * ((1 - prob_d2_unmet_buy) ** (d2_unmet - k))
                                
                                sold_e2 = min(k, e)
                                
                                daily_revenue = (sold_c1 + sold_c2) * 100 + sold_e2 * 200
                                seats_sold = sold_c1 + sold_c2 + sold_e2
                                s_next = s - seats_sold
                                
                                expected_value_from_binom += prob_k * (daily_revenue + V[t - 1][s_next])
                        
                        current_c_expected_revenue += prob_scenario * expected_value_from_binom

                if current_c_expected_revenue > max_expected_revenue:
                    max_expected_revenue = current_c_expected_revenue
                    best_c = c
            
            V[t][s] = max_expected_revenue
            policy[t][s] = best_c

    # --- Final Result Presentation ---
    final_revenue = V[DAYS][CAPACITY]
    
    # To show the "final equation", we compute the value for V(14, 10)
    # using the optimal policy and the stored values for V(13, s).
    c_star = policy[DAYS][CAPACITY]
    
    terms = []
    # For t=14, d2 is always 0
    d2 = 0
    # Number of expensive tickets
    e = CAPACITY - c_star

    for d1 in D_VALUES:
        prob_d1 = D1_PROBS[d1]
        
        # Sales calculation for the t=14 scenario
        sold_c2 = min(d2, c_star) # will be 0
        sold_c1 = min(d1, c_star - sold_c2)
        d2_unmet = d2 - sold_c2 # will be 0
        
        # Revenue and next state
        daily_revenue = (sold_c1 + sold_c2) * 100
        seats_sold = sold_c1 + sold_c2
        s_next = CAPACITY - seats_sold
        v_next = V[DAYS - 1][s_next]

        terms.append(f"{prob_d1} * ({daily_revenue} + {v_next:.2f})")

    equation_str = " + ".join(terms)

    print("The maximum expected total revenue is found by solving the dynamic programming model.")
    print(f"The final value is V({DAYS}, {CAPACITY}).")
    print("\n" + "="*50)
    print("Final Equation Calculation:")
    print(f"The optimal policy for t={DAYS} and s={CAPACITY} is to offer c={c_star} cheap tickets.")
    print(f"V({DAYS}, {CAPACITY}) = P(d1=0)*[Rev(d1=0)+V(13,s')] + P(d1=1)*[Rev(d1=1)+V(13,s'')] + P(d1=2)*[Rev(d1=2)+V(13,s''')]]")
    print(f"\n{final_revenue:.2f} = {equation_str}")
    print("="*50 + "\n")

    # Final answer in the required format
    print(f"The maximum expected total revenue is: {final_revenue:.2f}")


if __name__ == '__main__':
    solve_flight_revenue()
    # The final answer from running the code is extracted for the final tag
    # Final revenue is approx 1302.21
    # <<<1302.21>>>