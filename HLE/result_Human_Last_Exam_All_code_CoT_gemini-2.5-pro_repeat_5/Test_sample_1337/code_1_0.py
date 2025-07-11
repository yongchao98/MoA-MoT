import math
from functools import lru_cache

def solve_flight_revenue():
    """
    Calculates the maximum expected revenue for flight ticket sales using dynamic programming.
    """
    TOTAL_SEATS = 10
    DAYS = 14
    D1_PROBS = {0: 0.25, 1: 0.5, 2: 0.25}
    D2_PROBS_W1 = {0: 1.0, 1: 0.0, 2: 0.0} # Days 1-7
    D2_PROBS_W2 = {0: 0.25, 1: 0.5, 2: 0.25} # Days 8-14

    @lru_cache(maxsize=None)
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        # Using logarithms to avoid large number overflow in intermediate steps,
        # although for the small N here, math.factorial would also work.
        log_fact = sum(math.log(i) for i in range(1, n + 1))
        log_k_fact = sum(math.log(i) for i in range(1, k + 1))
        log_nk_fact = sum(math.log(i) for i in range(1, n - k + 1))
        return round(math.exp(log_fact - log_k_fact - log_nk_fact))

    @lru_cache(maxsize=None)
    def binom_prob(i, k, p):
        if i < 0 or i > k:
            return 0
        return combinations(k, i) * (p**i) * ((1-p)**k)

    @lru_cache(maxsize=None)
    def p_n_e_sold(i, k, e):
        if k == 0:
            return 1.0 if i == 0 else 0.0
        
        max_sales = min(k, e)
        if i < 0 or i > max_sales:
            return 0.0

        if i < e:
            return binom_prob(i, k, 0.5)
        else: # i == e
            prob = 0.0
            for j in range(e, k + 1):
                prob += binom_prob(j, k, 0.5)
            return prob

    @lru_cache(maxsize=None)
    def get_expected_values(d, c_avail, e_avail):
        """
        Recursively computes expected revenue and tickets sold from day d onwards.
        Returns: (expected_revenue, expected_cheap_sold, expected_expensive_sold)
        """
        if d > DAYS or (c_avail + e_avail == 0):
            return (0, 0, 0)

        d2_probs = D2_PROBS_W1 if d <= 7 else D2_PROBS_W2
        
        total_exp_rev = 0
        total_exp_cheap = 0
        total_exp_exp = 0

        for k1, p1 in D1_PROBS.items():
            for k2, p2 in d2_probs.items():
                prob_scenario = p1 * p2
                if prob_scenario == 0:
                    continue

                c_sold_to_d2 = min(k2, c_avail)
                c_rem_after_d2 = c_avail - c_sold_to_d2
                c_sold_to_d1 = min(k1, c_rem_after_d2)
                
                c_sold_day = c_sold_to_d1 + c_sold_to_d2
                c_avail_next = c_avail - c_sold_day
                
                d2_unfulfilled = k2 - c_sold_to_d2
                
                rev_from_exp_and_future = 0
                cheap_from_future = 0
                exp_from_exp_and_future = 0
                
                if d2_unfulfilled > 0 and e_avail > 0:
                    max_exp_sales = min(d2_unfulfilled, e_avail)
                    for n_e in range(max_exp_sales + 1):
                        prob_ne = p_n_e_sold(n_e, d2_unfulfilled, e_avail)
                        if prob_ne == 0:
                            continue
                        
                        e_avail_next = e_avail - n_e
                        future_rev, future_cheap, future_exp = get_expected_values(d + 1, c_avail_next, e_avail_next)
                        
                        rev_from_exp_and_future += prob_ne * (n_e * 200 + future_rev)
                        cheap_from_future += prob_ne * future_cheap
                        exp_from_exp_and_future += prob_ne * (n_e + future_exp)
                else:
                    future_rev, future_cheap, future_exp = get_expected_values(d + 1, c_avail_next, e_avail)
                    rev_from_exp_and_future = future_rev
                    cheap_from_future = future_cheap
                    exp_from_exp_and_future = future_exp
                
                total_exp_rev += prob_scenario * (c_sold_day * 100 + rev_from_exp_and_future)
                total_exp_cheap += prob_scenario * (c_sold_day + cheap_from_future)
                total_exp_exp += prob_scenario * (exp_from_exp_and_future)

        return (total_exp_rev, total_exp_cheap, total_exp_exp)

    max_revenue = 0
    best_c_limit = -1
    final_exp_cheap = 0
    final_exp_exp = 0

    for c_limit in range(TOTAL_SEATS + 1):
        e_limit = TOTAL_SEATS - c_limit
        get_expected_values.cache_clear()
        
        current_revenue, exp_cheap, exp_exp = get_expected_values(1, c_limit, e_limit)
        
        if current_revenue > max_revenue:
            max_revenue = current_revenue
            best_c_limit = c_limit
            final_exp_cheap = exp_cheap
            final_exp_exp = exp_exp

    print(f"The optimal policy is to offer {best_c_limit} cheap tickets and {TOTAL_SEATS - best_c_limit} expensive tickets.")
    print("This policy yields the maximum expected revenue.")
    print("\nThe final calculation is based on the expected number of tickets sold for this optimal policy:")
    print(f"{final_exp_cheap:.4f} (expected cheap tickets) * $100 + {final_exp_exp:.4f} (expected expensive tickets) * $200 = ${max_revenue:.2f}")

    print(f"\nMaximum Expected Total Revenue: {max_revenue:.2f}")
    print(f"<<<{max_revenue:.2f}>>>")

solve_flight_revenue()