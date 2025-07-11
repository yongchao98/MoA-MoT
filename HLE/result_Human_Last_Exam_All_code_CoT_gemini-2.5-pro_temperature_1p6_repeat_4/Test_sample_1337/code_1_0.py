import math
from functools import lru_cache

def solve_flight_revenue():
    """
    Calculates the maximum expected revenue for flight ticket sales
    using dynamic programming.
    """

    # --- Pre-computation for binomial probabilities ---
    # Used for Class 2 customers deciding to buy expensive tickets.
    memo_combinations = {}
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        if (n, k) in memo_combinations:
            return memo_combinations[(n, k)]
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        memo_combinations[(n, k)] = res
        return res

    memo_binom_pmf = {}
    def binom_pmf(k, n, p):
        # Probability Mass Function for a Binomial distribution
        if (k, n, p) in memo_binom_pmf:
            return memo_binom_pmf[(k, n, p)]
        if k < 0 or k > n:
            return 0.0
        res = combinations(n, k) * (p**k) * ((1-p)**(n-k))
        memo_binom_pmf[(k, n, p)] = res
        return res

    # --- Demand Distributions as per the problem description ---
    P_D1 = {0: 0.25, 1: 0.5, 2: 0.25}  # For Class 1 customers (all 14 days)
    P_D2_WEEK1 = {0: 1.0}              # For Class 2 customers in week 1 (days 1-7)
    P_D2_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25} # For Class 2 customers in week 2 (days 8-14)
    
    # --- Main DP function with memoization ---
    @lru_cache(maxsize=None)
    def calculate_expected_values(t, c_rem, s_rem):
        """
        Calculates expected outcomes from a given state.
        Args:
            t (int): The current day (1 to 14).
            c_rem (int): The remaining quota of cheap tickets.
            s_rem (int): The remaining quota of expensive tickets.
        Returns:
            tuple: (expected_future_revenue, exp_future_cheap_sales, exp_future_expensive_sales)
        """
        # Base case: After day 14 or if no seats are left, future revenue is zero.
        if t > 14 or (c_rem + s_rem == 0):
            return (0.0, 0.0, 0.0)

        # Determine Class 2 demand distribution for the current week
        P_D2 = P_D2_WEEK1 if 1 <= t <= 7 else P_D2_WEEK2

        # Initialize total expected values for the current state
        total_exp_rev = 0.0
        total_exp_c_sold = 0.0
        total_exp_s_sold = 0.0

        # Iterate over all possible demand scenarios for the day
        for d1, p_d1 in P_D1.items():
            for d2, p_d2 in P_D2.items():
                prob_demand = p_d1 * p_d2
                if prob_demand == 0:
                    continue
                
                # --- Simulate ticket sales for this demand scenario (d1, d2) ---
                # Class 2 customers have priority for cheap tickets
                c2_buy_cheap = min(d2, c_rem)
                c_after_c2 = c_rem - c2_buy_cheap
                # Class 1 customers take the rest
                c1_buy_cheap = min(d1, c_after_c2)
                
                c_sold_today = c1_buy_cheap + c2_buy_cheap
                rev_cheap_today = c_sold_today * 100
                c_new = c_rem - c_sold_today

                # Class 2 customers who couldn't get a cheap ticket might buy an expensive one
                d2_unhappy = d2 - c2_buy_cheap
                exp_seats_avail = s_rem
                
                # These will be the expected values calculated for this specific demand (d1, d2)
                exp_rev_demand, exp_c_demand, exp_s_demand = (0.0, 0.0, 0.0)

                # Iterate over possible number of expensive sales `k`
                # This follows a Binomial distribution B(d2_unhappy, 0.5), capped by exp_seats_avail.
                for k in range(min(d2_unhappy, exp_seats_avail) + 1):
                    # Calculate the probability of k expensive sales
                    prob_k_sales = 0.0
                    if k < exp_seats_avail:
                        prob_k_sales = binom_pmf(k, d2_unhappy, 0.5)
                    else: # k == exp_seats_avail, we sold all remaining expensive seats
                        for j in range(k, d2_unhappy + 1):
                            prob_k_sales += binom_pmf(j, d2_unhappy, 0.5)

                    if prob_k_sales == 0:
                        continue

                    rev_exp_today = k * 200
                    s_new = s_rem - k
                    
                    # Recursively get expected values from the next day onwards
                    future_rev, future_c, future_s = calculate_expected_values(t + 1, c_new, s_new)
                    
                    # Add this outcome's contribution to the expected values for this demand
                    exp_rev_demand += prob_k_sales * (rev_cheap_today + rev_exp_today + future_rev)
                    exp_c_demand += prob_k_sales * (c_sold_today + future_c)
                    exp_s_demand += prob_k_sales * (k + future_s)
                
                # Add the weighted results for this demand to the daily total
                total_exp_rev += prob_demand * exp_rev_demand
                total_exp_c_sold += prob_demand * exp_c_demand
                total_exp_s_sold += prob_demand * exp_s_demand

        return (total_exp_rev, total_exp_c_sold, total_exp_s_sold)

    # --- Main loop to find the optimal policy ---
    max_rev = 0.0
    optimal_c_limit = -1
    optimal_sales_breakdown = (0, 0)
    
    # Iterate through all possible static policies (number of cheap tickets to allocate)
    for c_limit in range(11):  # c_limit can be 0 to 10
        s_limit = 10 - c_limit
        calculate_expected_values.cache_clear()
        
        # Calculate expected outcomes for this policy
        current_rev, exp_c, exp_s = calculate_expected_values(1, c_limit, s_limit)

        if current_rev > max_rev:
            max_rev = current_rev
            optimal_c_limit = c_limit
            optimal_sales_breakdown = (exp_c, exp_s)

    # --- Final Output ---
    opt_exp_c_sold, opt_exp_s_sold = optimal_sales_breakdown
    
    print(f"The optimal policy is to allocate {optimal_c_limit} tickets as cheap and {10-optimal_c_limit} as expensive.")
    print("This yields the maximum expected total revenue.")
    print("\nFinal equation based on the optimal policy:")
    print(f"E[Revenue] = E[Cheap Sales] * $100 + E[Expensive Sales] * $200")
    print(f"{max_rev:.2f} = {opt_exp_c_sold:.4f} * 100 + {opt_exp_s_sold:.4f} * 200")

    # The final answer in the required format
    print(f"\nMaximum Expected Total Revenue: ${max_rev:.2f}")
    print(f"<<<{max_rev:.2f}>>>")

solve_flight_revenue()