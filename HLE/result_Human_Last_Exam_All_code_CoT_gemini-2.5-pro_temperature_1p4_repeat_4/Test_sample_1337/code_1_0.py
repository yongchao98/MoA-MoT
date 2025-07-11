import math
from functools import lru_cache

def solve_revenue_management():
    """
    Solves the airline revenue management problem using dynamic programming.
    """
    # --- Problem Constants ---
    SEAT_CAPACITY = 10
    DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200

    # --- Customer Arrival Probabilities ---
    P1_dist = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1
    P2_dist = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2

    # Memoization for binomial probability: P(k successes in N trials with probability p)
    @lru_cache(maxsize=None)
    def binomial_prob(k, N, p=0.5):
        if k < 0 or k > N:
            return 0
        return math.comb(N, k) * (p ** k) * ((1 - p) ** (N - k))

    # --- Main Loop to Find the Best Policy (c_limit) ---
    max_revenue = -1.0
    best_policy_details = None

    for c_limit in range(SEAT_CAPACITY + 1):
        # DP table: dp[t][s_c][s_e] stores a tuple:
        # (expected_revenue, expected_cheap_sales, expected_expensive_sales)
        # Dimensions: t=days+2, s_c=seats+1, s_e=seats+1 for simplicity
        dp = [[[(0.0, 0.0, 0.0) for _ in range(SEAT_CAPACITY + 1)] 
               for _ in range(SEAT_CAPACITY + 1)] 
              for _ in range(DAYS + 2)]

        # --- Dynamic Programming: Backward Induction ---
        for t in range(DAYS, 0, -1):  # From day 14 down to day 1
            for s_c in range(c_limit + 1):
                for s_e in range(SEAT_CAPACITY - s_c + 1):
                    
                    if s_c + s_e >= SEAT_CAPACITY:
                        continue  # No more seats to sell

                    # Expected values for this state (t, s_c, s_e)
                    state_exp_rev, state_exp_cheap, state_exp_exp = 0.0, 0.0, 0.0

                    # Iterate over all possible customer arrival scenarios for day t
                    for i in P1_dist:  # Class 1 arrivals
                        for j in P2_dist:  # Class 2 arrivals
                            
                            if t <= 7: # First week: no Class 2 customers
                                if j > 0: continue
                                prob_arrival = P1_dist[i]
                            else: # Second week: both classes
                                prob_arrival = P1_dist[i] * P2_dist[j]
                            
                            # --- Sales Logic for one arrival scenario (i, j) ---
                            seats_left = SEAT_CAPACITY - (s_c + s_e)
                            cheap_tickets_policy_left = c_limit - s_c
                            
                            cheap_available_for_sale = min(seats_left, cheap_tickets_policy_left)
                            
                            c2_gets_cheap = min(j, cheap_available_for_sale)
                            c1_gets_cheap = min(i, cheap_available_for_sale - c2_gets_cheap)
                            sold_c_today = c1_gets_cheap + c2_gets_cheap

                            # --- Expensive Ticket Sales Logic ---
                            c2_no_cheap = j - c2_gets_cheap
                            expensive_seats_available = seats_left - sold_c_today
                            
                            # Expected values, conditioned on arrival (i,j)
                            # We take an expectation over the choices of the 'c2_no_cheap' customers
                            exp_val_cond_on_arrival = (0.0, 0.0, 0.0)
                            
                            N = c2_no_cheap  # Num customers considering expensive ticket
                            M = expensive_seats_available # Num expensive seats available

                            for k in range(N + 1): # k is the number who decide to buy
                                prob_k_buys = binomial_prob(k, N)
                                sold_e_today = min(k, M)
                                
                                revenue_today = sold_c_today * PRICE_CHEAP + sold_e_today * PRICE_EXPENSIVE
                                
                                # Get expected future values from dp table
                                future_rev, future_cheap, future_exp = dp[t + 1][s_c + sold_c_today][s_e + sold_e_today]
                                
                                # Summing the expectations
                                r, c, e = exp_val_cond_on_arrival
                                r += prob_k_buys * (revenue_today + future_rev)
                                c += prob_k_buys * (sold_c_today + future_cheap)
                                e += prob_k_buys * (sold_e_today + future_exp)
                                exp_val_cond_on_arrival = (r, c, e)

                            # Update total expected value for the state
                            state_exp_rev += prob_arrival * exp_val_cond_on_arrival[0]
                            state_exp_cheap += prob_arrival * exp_val_cond_on_arrival[1]
                            state_exp_exp += prob_arrival * exp_val_cond_on_arrival[2]

                    dp[t][s_c][s_e] = (state_exp_rev, state_exp_cheap, state_exp_exp)
        
        # Result for the current c_limit policy
        initial_exp_rev, initial_exp_cheap, initial_exp_exp = dp[1][0][0]
        
        if initial_exp_rev > max_revenue:
            max_revenue = initial_exp_rev
            best_policy_details = (initial_exp_rev, initial_exp_cheap, initial_exp_exp)

    # --- Print the Final Result ---
    if best_policy_details:
        final_rev, final_cheap_sales, final_expensive_sales = best_policy_details
        print(f"{final_cheap_sales:.4f} * {PRICE_CHEAP} + {final_expensive_sales:.4f} * {PRICE_EXPENSIVE} = {final_rev:.4f}")

# Execute the solver
solve_revenue_management()