import math

def solve_flight_revenue():
    """
    Solves the flight revenue management problem using dynamic programming.
    """
    # --- Problem Parameters ---
    CAPACITY = 10
    DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200

    # Demand Distributions P(k requests)
    P1 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1
    P2 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2

    # --- DP Initialization ---
    # V[d][s]: max expected revenue with d days left and s seats
    V = [[0.0 for _ in range(CAPACITY + 1)] for _ in range(DAYS + 1)]
    # policy[d][s]: optimal protection level p for state (d, s)
    policy = [[0 for _ in range(CAPACITY + 1)] for _ in range(DAYS + 1)]

    # Memoization for binomial probabilities of C2 expensive purchase (p=0.5)
    binom_memo = {}
    def p_binom(k, n):
        if (k, n) in binom_memo:
            return binom_memo[(k, n)]
        if k < 0 or k > n:
            return 0
        prob = math.comb(n, k) * (0.5 ** n)
        binom_memo[(k, n)] = prob
        return prob

    # --- DP Calculation ---
    # d = days remaining
    for d in range(1, DAYS + 1):
        # s = seats available
        for s in range(CAPACITY + 1):
            max_ev_for_s = -1.0
            optimal_p_for_s = 0

            # p = protection level (seats reserved for expensive fare)
            for p in range(s + 1):
                c_avail_total = s - p
                current_p_ev = 0.0

                # Week 1: days 14 down to 8 (d > 7). No Class 2 customers.
                if d > 7:
                    for d1, prob_d1 in P1.items():
                        c1_sold = min(d1, c_avail_total)
                        revenue = c1_sold * PRICE_CHEAP
                        s_next = s - c1_sold
                        current_p_ev += prob_d1 * (revenue + V[d - 1][s_next])
                
                # Week 2: days 7 down to 1 (d <= 7). Both customer classes.
                else:
                    for d1, prob_d1 in P1.items():
                        for d2, prob_d2 in P2.items():
                            prob_path = prob_d1 * prob_d2

                            # Sales of cheap tickets (C2 has priority)
                            c2_sold = min(d2, c_avail_total)
                            c_avail_for_c1 = c_avail_total - c2_sold
                            c1_sold = min(d1, c_avail_for_c1)
                            rev_cheap = (c1_sold + c2_sold) * PRICE_CHEAP
                            
                            s_after_cheap = s - c1_sold - c2_sold
                            d2_rem = d2 - c2_sold # C2 customers who might buy expensive

                            # Expected value from expensive sales stage onwards
                            ev_from_exp_stage = 0.0
                            if d2_rem == 0:
                                ev_from_exp_stage = V[d - 1][s_after_cheap]
                            else:
                                for k_buy in range(d2_rem + 1):
                                    prob_k_buy = p_binom(k_buy, d2_rem)
                                    e2_sold = min(k_buy, s_after_cheap)
                                    rev_exp = e2_sold * PRICE_EXPENSIVE
                                    s_next = s_after_cheap - e2_sold
                                    ev_from_exp_stage += prob_k_buy * (rev_exp + V[d - 1][s_next])
                            
                            current_p_ev += prob_path * (rev_cheap + ev_from_exp_stage)

                # Update max EV and policy for state (d, s)
                if current_p_ev > max_ev_for_s:
                    max_ev_for_s = current_p_ev
                    optimal_p_for_s = p
            
            V[d][s] = max_ev_for_s
            policy[d][s] = optimal_p_for_s

    # --- Final Result ---
    max_revenue = V[DAYS][CAPACITY]
    p_star = policy[DAYS][CAPACITY]
    c_avail_star = CAPACITY - p_star

    # Calculate terms for the "final equation" for the first day's decision
    terms = []
    total_check = 0
    # Day 14 is in the first week, so only Class 1 demand matters for the first step
    for d1, prob_d1 in P1.items():
        c1_sold = min(d1, c_avail_star)
        revenue = c1_sold * PRICE_CHEAP
        s_next = CAPACITY - c1_sold
        future_v = V[DAYS - 1][s_next]
        term_value = prob_d1 * (revenue + future_v)
        total_check += term_value
        terms.append(f"{prob_d1} * ({revenue} + {future_v:.2f})")

    print("The maximum expected total revenue is determined by solving a dynamic programming model.")
    print("Let V(d, s) be the max expected revenue with d days left and s seats.")
    print("The optimal policy is to choose a protection level 'p' to maximize V(d, s).")
    print(f"\nFor the initial state (14 days, 10 seats), the optimal protection level is p* = {p_star}")
    print(f"This means we should initially offer {c_avail_star} cheap tickets and protect {p_star} seats for high-fare customers.")
    print("\nThe calculation for V(14, 10), based on the first day's events, is:")
    equation = f"V(14, 10) = " + " + ".join(terms)
    print(equation)
    print(f"            = {total_check:.2f}")

    print(f"\nThe maximum expected total revenue is: {max_revenue:.2f}")
    
    return max_revenue

if __name__ == '__main__':
    solve_flight_revenue()