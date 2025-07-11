import math
from functools import lru_cache

def solve_revenue_management():
    """
    Calculates the maximum expected revenue using dynamic programming.
    """
    # --- Problem Parameters ---
    CAPACITY = 10
    DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200
    DEMAND_PROBS = {0: 0.25, 1: 0.5, 2: 0.25}

    # DP table V[t][s] stores the max expected revenue with t days left and s seats.
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(DAYS + 1)]

    # --- Main DP Loop ---
    # Iterate backwards in time, from 1 day left to 14 days left.
    for t in range(1, DAYS + 1):

        # This memoized helper function calculates the expected value from
        # the probabilistic sale of expensive tickets.
        # It's redefined inside the 't' loop to correctly capture V[t-1].
        @lru_cache(maxsize=None)
        def get_exp_value_from_exp_sales(buyers, seats):
            # Base case: no more potential buyers or no seats left.
            # The only remaining value is the future value from the next state.
            if buyers == 0 or seats == 0:
                return V[t - 1][seats]

            # Recursive step: consider the first potential buyer.
            # Case 1: Buyer purchases with 50% probability.
            # Revenue increases by the expensive price, and one buyer and one seat are used.
            val_if_buy = PRICE_EXPENSIVE + get_exp_value_from_exp_sales(buyers - 1, seats - 1)

            # Case 2: Buyer walks away with 50% probability.
            # No revenue is gained, one buyer is used, but the number of seats remains.
            val_if_walk = 0 + get_exp_value_from_exp_sales(buyers - 1, seats)
            
            # The expected value is the average of the two outcomes.
            return 0.5 * val_if_buy + 0.5 * val_if_walk

        # Iterate over all possible numbers of available seats 's'.
        for s in range(1, CAPACITY + 1):
            max_rev_for_ts = -1.0
            
            # Iterate over all possible protection levels 'p'.
            for p in range(s + 1):
                current_p_expected_rev = 0.0
                
                # Determine demand distribution for Class 2 customers based on the week.
                is_second_week = (t <= 7)
                demand2_dist = DEMAND_PROBS if is_second_week else {0: 1.0}
                
                # Sum expectation over all demand scenarios for Class 1 (d1) and Class 2 (d2).
                for d1, p1 in DEMAND_PROBS.items():
                    for d2, p2 in demand2_dist.items():
                        prob_scenario = p1 * p2
                        
                        # --- Simulate one day's sales for this policy and demand ---
                        seats_cheap_avail = max(0, s - p)
                        
                        # 1. Class 2 customers buy cheap tickets first (priority).
                        sold2_cheap = min(d2, seats_cheap_avail)
                        rev_c2_cheap = sold2_cheap * PRICE_CHEAP
                        seats_rem_after_c2c = s - sold2_cheap
                        seats_cheap_rem_after_c2c = seats_cheap_avail - sold2_cheap
                        d2_unmet_cheap = d2 - sold2_cheap
                        
                        # 2. Class 1 customers buy remaining cheap tickets.
                        sold1_cheap = min(d1, seats_cheap_rem_after_c2c)
                        rev_c1_cheap = sold1_cheap * PRICE_CHEAP
                        seats_rem_after_c1 = seats_rem_after_c2c - sold1_cheap
                        
                        day_rev_so_far = rev_c1_cheap + rev_c2_cheap
                        
                        # 3. Unmet Class 2 demand considers expensive tickets.
                        potential_exp_buyers = d2_unmet_cheap
                        seats_for_exp = seats_rem_after_c1
                        
                        # This call returns E[day_rev_from_exp_sales + future_rev]
                        exp_val_from_exp_part = get_exp_value_from_exp_sales(
                            potential_exp_buyers, seats_for_exp
                        )
                        
                        scenario_total_value = day_rev_so_far + exp_val_from_exp_part
                        current_p_expected_rev += prob_scenario * scenario_total_value

                # Update the max revenue for state (t,s) if this policy 'p' is better.
                if current_p_expected_rev > max_rev_for_ts:
                    max_rev_for_ts = current_p_expected_rev
            
            V[t][s] = max_rev_for_ts
        
        # Clear the cache for the helper function before the next time step.
        get_exp_value_from_exp_sales.cache_clear()

    # The final answer is in V[14][10].
    final_result = V[DAYS][CAPACITY]
    
    # Print the result in a descriptive way.
    print(f"The calculation for a flight with a capacity of {CAPACITY} seats over a {DAYS} day selling period is complete.")
    print(f"The maximum expected total revenue is: ${final_result:.2f}")

solve_revenue_management()
<<<1357.94>>>