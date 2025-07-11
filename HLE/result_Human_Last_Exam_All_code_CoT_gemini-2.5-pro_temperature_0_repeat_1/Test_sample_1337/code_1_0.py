import collections

def solve():
    """
    Calculates the maximum expected total revenue using dynamic programming.
    """
    CAPACITY = 10
    DAYS = 14

    # Customer demand distributions
    # P1: Class 1 customers (Days 1-14)
    P1 = {0: 0.25, 1: 0.5, 2: 0.25}
    # P2_WEEK1: Class 2 customers (Days 1-7)
    P2_WEEK1 = {0: 1.0, 1: 0.0, 2: 0.0}
    # P2_WEEK2: Class 2 customers (Days 8-14)
    P2_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25}

    # V[t, s]: max expected revenue from day t onwards with s seats left.
    # Using a defaultdict simplifies handling the base case (V[15, s] = 0).
    V = collections.defaultdict(float)

    # Loop backwards from the last day of sales to the first
    for t in range(DAYS, 0, -1):
        # Determine the correct demand distribution for Class 2 customers
        P2 = P2_WEEK1 if t <= 7 else P2_WEEK2
        
        # Loop through all possible numbers of remaining seats
        for s in range(CAPACITY + 1):
            max_rev_for_state = 0.0
            
            # The decision: how many of the s seats to offer as cheap (c_policy)
            for c_policy in range(s + 1):
                e_policy = s - c_policy
                
                expected_rev_for_this_policy = 0.0
                # Iterate over all possible demand scenarios for this day
                for d1, p1 in P1.items():
                    for d2, p2 in P2.items():
                        prob = p1 * p2
                        if prob == 0:
                            continue

                        # --- Simulate sales for the day ---
                        immediate_rev = 0.0
                        
                        # 1. Class 2 customers have priority for cheap tickets
                        sold_c2_cheap = min(d2, c_policy)
                        immediate_rev += sold_c2_cheap * 100
                        d2_rem = d2 - sold_c2_cheap

                        # 2. Class 1 customers buy remaining cheap tickets
                        c_rem_for_c1 = c_policy - sold_c2_cheap
                        sold_c1_cheap = min(d1, c_rem_for_c1)
                        immediate_rev += sold_c1_cheap * 100
                        
                        seats_sold_cheap = sold_c1_cheap + sold_c2_cheap
                        
                        # 3. Remaining Class 2 customers consider expensive tickets.
                        # We calculate the expected value, which is the sum of immediate revenue
                        # from expensive sales and the expected future revenue V(t+1, s').
                        exp_c2_exp_and_future_rev = 0.0
                        
                        if d2_rem == 0:
                            s_prime = s - seats_sold_cheap
                            exp_c2_exp_and_future_rev = V[t + 1, s_prime]
                        elif d2_rem == 1:
                            if e_policy >= 1:
                                # Buys with 50% prob, walks away with 50% prob
                                s_prime_sale = s - seats_sold_cheap - 1
                                s_prime_no_sale = s - seats_sold_cheap
                                rev_if_sale = 200 + V[t + 1, s_prime_sale]
                                rev_if_no_sale = V[t + 1, s_prime_no_sale]
                                exp_c2_exp_and_future_rev = 0.5 * rev_if_sale + 0.5 * rev_if_no_sale
                            else: # No expensive seats available
                                s_prime = s - seats_sold_cheap
                                exp_c2_exp_and_future_rev = V[t + 1, s_prime]
                        elif d2_rem == 2:
                            if e_policy == 0:
                                s_prime = s - seats_sold_cheap
                                exp_c2_exp_and_future_rev = V[t + 1, s_prime]
                            elif e_policy == 1:
                                # 1 sale with prob 0.75, 0 sales with prob 0.25
                                s_prime_1_sale = s - seats_sold_cheap - 1
                                s_prime_0_sales = s - seats_sold_cheap
                                rev_if_1_sale = 200 + V[t + 1, s_prime_1_sale]
                                rev_if_0_sales = V[t + 1, s_prime_0_sales]
                                exp_c2_exp_and_future_rev = 0.75 * rev_if_1_sale + 0.25 * rev_if_0_sales
                            else:  # e_policy >= 2
                                # 2 sales (p=0.25), 1 sale (p=0.5), 0 sales (p=0.25)
                                s_prime_2_sales = s - seats_sold_cheap - 2
                                s_prime_1_sale = s - seats_sold_cheap - 1
                                s_prime_0_sales = s - seats_sold_cheap
                                rev_if_2_sales = 400 + V[t + 1, s_prime_2_sales]
                                rev_if_1_sale = 200 + V[t + 1, s_prime_1_sale]
                                rev_if_0_sales = V[t + 1, s_prime_0_sales]
                                exp_c2_exp_and_future_rev = 0.25 * rev_if_2_sales + 0.5 * rev_if_1_sale + 0.25 * rev_if_0_sales
                        
                        total_exp_rev_for_scenario = immediate_rev + exp_c2_exp_and_future_rev
                        expected_rev_for_this_policy += prob * total_exp_rev_for_scenario
                
                # Update the max revenue for the current state (t,s)
                if expected_rev_for_this_policy > max_rev_for_state:
                    max_rev_for_state = expected_rev_for_this_policy
            
            V[t, s] = max_rev_for_state

    # The final answer is the max expected revenue at the start (day 1, 10 seats)
    max_expected_revenue = V[1, CAPACITY]
    
    # The final equation is simply the maximum expected revenue value.
    print(f"The maximum expected total revenue is: {max_expected_revenue:.2f}")

solve()