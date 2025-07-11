import math

def solve():
    """
    Calculates the maximum expected revenue using dynamic programming.
    """
    CAPACITY = 10
    DAYS = 14

    # DP table V[t][s]: max expected revenue with t days left and s seats.
    # V has dimensions (DAYS + 1) x (CAPACITY + 1).
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(DAYS + 1)]

    # Pre-computation for Class 2 expensive ticket purchases.
    # E_rev[n][s]: expected revenue from n customers with s seats.
    # E_s_sold[n][s]: expected seats sold to n customers with s seats.
    # Maximum Class 2 demand is 2, so we need tables of size 3.
    E_rev = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(3)]
    E_s_sold = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(3)]

    # Base cases for n=0 (no customers) are already 0.
    # Compute for n=1 customer.
    for s in range(1, CAPACITY + 1):
        # A single customer buys an expensive ticket ($200) with 50% probability.
        E_rev[1][s] = 0.5 * 200.0
        E_s_sold[1][s] = 0.5 * 1.0

    # Compute for n=2 customers using the results for n=1.
    for s in range(1, CAPACITY + 1):
        # The recurrence for expected revenue:
        E_rev[2][s] = 0.5 * (200.0 + E_rev[1][s - 1]) + 0.5 * E_rev[1][s]
        # The recurrence for expected seats sold:
        E_s_sold[2][s] = 0.5 * (1.0 + E_s_sold[1][s - 1]) + 0.5 * E_s_sold[1][s]

    # Customer arrival probability distributions.
    P_D1 = {0: 0.25, 1: 0.5, 2: 0.25}
    P_D2_first_week = {0: 1.0, 1: 0.0, 2: 0.0} # Days 1-7
    P_D2_second_week = {0: 0.25, 1: 0.5, 2: 0.25} # Days 8-14

    # Main DP loop, working backward from departure.
    # t: days remaining before departure.
    for t in range(1, DAYS + 1):
        day_of_sale = DAYS - t + 1
        P_D2 = P_D2_first_week if day_of_sale <= 7 else P_D2_second_week

        # s: number of seats available.
        for s in range(1, CAPACITY + 1):
            max_rev_for_s = -1.0
            
            # k: protection level (number of seats reserved for expensive fares).
            for k in range(s + 1):
                # c: number of cheap tickets offered today.
                c = s - k 
                
                current_k_total_expected_revenue = 0.0
                
                # Iterate over all combinations of customer arrivals.
                for i, p_i in P_D1.items():
                    for j, p_j in P_D2.items():
                        prob = p_i * p_j
                        if prob == 0:
                            continue
                        
                        # --- Simulate sales for this (i, j) arrival outcome ---
                        
                        # Class 2 customers have priority for cheap tickets.
                        cheap_sold_c2 = min(j, c)
                        
                        # Class 1 customers try for remaining cheap tickets.
                        cheap_seats_left = c - cheap_sold_c2
                        cheap_sold_c1 = min(i, cheap_seats_left)
                        
                        rev_cheap = (cheap_sold_c1 + cheap_sold_c2) * 100.0
                        
                        # Remaining C2 customers who might buy expensive tickets.
                        j_rem = j - cheap_sold_c2
                        seats_after_cheap = s - cheap_sold_c1 - cheap_sold_c2
                        
                        # Expected revenue and sales from these C2 customers.
                        rev_exp = E_rev[j_rem][seats_after_cheap]
                        exp_sold_count = E_s_sold[j_rem][seats_after_cheap]
                        
                        day_revenue = rev_cheap + rev_exp
                        
                        # Expected number of seats left for the future.
                        s_final_expected = seats_after_cheap - exp_sold_count
                        
                        # --- Calculate expected future revenue using interpolation ---
                        s_low = math.floor(s_final_expected)
                        s_high = math.ceil(s_final_expected)
                        
                        future_rev = 0.0
                        if s_low >= 0:
                            if s_low == s_high:
                                future_rev = V[t - 1][int(s_low)]
                            else:
                                w_low = s_high - s_final_expected
                                w_high = s_final_expected - s_low
                                future_rev = V[t-1][int(s_low)] * w_low + V[t-1][int(s_high)] * w_high
                        
                        current_k_total_expected_revenue += prob * (day_revenue + future_rev)

                if current_k_total_expected_revenue > max_rev_for_s:
                    max_rev_for_s = current_k_total_expected_revenue

            V[t][s] = max_rev_for_s
            
    # The final answer is the max expected revenue with 10 seats and 14 days to go.
    final_revenue = V[DAYS][CAPACITY]
    
    print(final_revenue)

solve()