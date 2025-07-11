import math

def solve_flight_revenue():
    """
    Solves the airline revenue management problem using dynamic programming.
    """

    # --- Setup ---
    TOTAL_DAYS = 14
    CAPACITY = 10
    CHEAP_PRICE = 100.0
    EXPENSIVE_PRICE = 200.0

    # DP table: V[d][s] = max expected revenue with d days left and s seats
    # Dimensions: (TOTAL_DAYS + 1) x (CAPACITY + 1)
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(TOTAL_DAYS + 1)]

    # Customer arrival probabilities for each class
    # P(arrivals=0)=0.25, P(arrivals=1)=0.5, P(arrivals=2)=0.25
    arrival_prob = {0: 0.25, 1: 0.5, 2: 0.25}

    # --- Helper for Binomial Probability ---
    # Used for C2 customers deciding to buy an expensive ticket
    comb_memo = {}
    def combinations(n, k):
        if k < 0 or k > n: return 0
        if (n, k) in comb_memo: return comb_memo[(n, k)]
        if k == 0 or k == n: return 1
        if k > n // 2: k = n - k
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        comb_memo[(n, k)] = res
        return res

    binom_prob_memo = {}
    def binom_prob(n, k, p):
        if (n,k,p) in binom_prob_memo: return binom_prob_memo[(n, k, p)]
        res = combinations(n, k) * (p**k) * ((1 - p)**(n - k))
        binom_prob_memo[(n,k,p)] = res
        return res


    # --- Dynamic Programming Calculation ---

    # Iterate backwards from d=1 to d=14
    for d in range(1, TOTAL_DAYS + 1):
        is_week_2 = (d <= 7)
        for s in range(1, CAPACITY + 1):
            max_exp_rev_for_s = 0.0
            
            # Decision: k = number of cheap tickets to offer today
            for k in range(s + 1):
                current_exp_rev_for_k = 0.0
                
                # Iterate over all arrival scenarios for n1 (Class 1) and n2 (Class 2)
                for n1, p1 in arrival_prob.items():
                    n2_arrivals = arrival_prob if is_week_2 else {0: 1.0}
                    for n2, p2 in n2_arrivals.items():
                        
                        prob_scenario = p1 * p2
                        
                        # --- Sales Logic ---
                        # C2 customers have priority for cheap tickets
                        sold_c2 = min(n2, k)
                        # C1 customers get the rest of the cheap tickets
                        sold_c1 = min(n1, k - sold_c2)
                        
                        revenue_from_cheap = CHEAP_PRICE * (sold_c1 + sold_c2)
                        
                        c2_turned_away = n2 - sold_c2
                        seats_left_for_exp = s - (sold_c1 + sold_c2)
                        
                        exp_val_from_expensive = 0.0
                        
                        # C2 customers turned away from cheap may buy expensive
                        if c2_turned_away > 0 and seats_left_for_exp > 0:
                            # Iterate over the number of C2 customers (m) who try to buy
                            for m in range(c2_turned_away + 1):
                                prob_m_buy = binom_prob(c2_turned_away, m, 0.5)
                                
                                sold_e = min(m, seats_left_for_exp)
                                revenue_from_expensive = EXPENSIVE_PRICE * sold_e
                                
                                s_next = s - (sold_c1 + sold_c2 + sold_e)
                                exp_val_from_expensive += prob_m_buy * (revenue_from_expensive + V[d - 1][s_next])
                        else:
                            # No expensive sales this turn
                            s_next = s - (sold_c1 + sold_c2)
                            exp_val_from_expensive = V[d-1][s_next]

                        current_exp_rev_for_k += prob_scenario * (revenue_from_cheap + exp_val_from_expensive)

                # Update the max expected revenue for state (d,s)
                if current_exp_rev_for_k > max_exp_rev_for_s:
                    max_exp_rev_for_s = current_exp_rev_for_k
            
            V[d][s] = max_exp_rev_for_s

    # --- Present the Final Result ---
    
    # Redo calculation for V(14, 10) to find the optimal k and show the breakdown
    d = 14
    s = 10
    opt_k = -1
    max_rev = -1.0
    
    for k in range(s + 1):
        exp_rev_for_k = 0.0
        for n1, p1 in arrival_prob.items():
            sold = min(n1, k)
            rev = CHEAP_PRICE * sold
            s_next = s - sold
            exp_rev_for_k += p1 * (rev + V[d-1][s_next])
        if exp_rev_for_k > max_rev:
            max_rev = exp_rev_for_k
            opt_k = k
    
    print("This problem is solved using dynamic programming. The maximum expected total revenue is V(14, 10).")
    print("\n--- Final Equation Breakdown ---")
    print(f"On the first day (d=14, s=10), the optimal number of cheap tickets to offer is k = {opt_k}.")
    print("The expected revenue is calculated based on the number of Class 1 arrivals (n1) for that day:")
    print("V(14, 10) = Σ [ P(n1) * (Revenue_today + V(13, s_next) ) ]")

    total_calc = 0
    full_eq_str = []
    
    for n1, p1 in arrival_prob.items():
        sold = min(n1, opt_k)
        rev = CHEAP_PRICE * sold
        s_next = s - sold
        v_future = V[d - 1][s_next]
        term = rev + v_future
        total_calc += p1 * term
        
        print(f"\nIf n1 = {n1} (Prob={p1}):")
        print(f"  Tickets sold: {sold}")
        print(f"  Revenue today: {rev}")
        print(f"  Seats remaining: {s_next}")
        print(f"  Future value V(13, {s_next}): {v_future:.2f}")
        print(f"  Term value: {p1} * ({rev:.0f} + {v_future:.2f}) = {p1 * term:.2f}")
        full_eq_str.append(f"{p1}*({rev:.0f}+{v_future:.2f})")

    print("\n--- Summary Calculation ---")
    print(" + ".join(full_eq_str) + f" = {total_calc:.4f}")

    final_revenue = V[TOTAL_DAYS][CAPACITY]
    print(f"\nMaximum Expected Total Revenue: {final_revenue:.4f}")
    
    # Print the final answer in the specified format
    print(f"<<<{final_revenue:.4f}>>>")

solve_flight_revenue()
