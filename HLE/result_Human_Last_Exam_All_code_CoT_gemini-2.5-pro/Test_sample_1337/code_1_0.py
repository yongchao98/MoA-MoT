import math
from functools import lru_cache

def solve_revenue_maximization():
    """
    This function sets up and solves the flight revenue maximization problem.
    """
    # --- Problem Constants ---
    CAPACITY = 10
    DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200
    PROB_C2_BUY_EXPENSIVE = 0.5

    # --- Customer Arrival Distributions ---
    DIST_C1 = {0: 0.25, 1: 0.5, 2: 0.25}
    DIST_C2_WEEK1 = {0: 1.0}
    DIST_C2_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25}

    # --- Helper function for binomial probability ---
    @lru_cache(maxsize=None)
    def binom_pmf(k, n, p):
        """
        Calculates the probability mass function of a binomial distribution.
        """
        if k < 0 or k > n:
            return 0
        return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

    # --- Main DP function with memoization ---
    @lru_cache(maxsize=None)
    def solve(t, c, e):
        """
        Calculates the maximum expected future revenue using dynamic programming.
        t: current day (from 1 to 14)
        c: number of cheap tickets available
        e: number of expensive tickets available
        """
        # Base case: selling period is over or no tickets are left
        if t > DAYS or (c == 0 and e == 0):
            return 0

        total_expected_value = 0
        
        # Determine Class 2 customer distribution based on the week
        dist_c2 = DIST_C2_WEEK1 if t <= 7 else DIST_C2_WEEK2

        # Iterate over all possible customer arrival scenarios for the day
        for d1, p1 in DIST_C1.items():
            for d2, p2 in dist_c2.items():
                prob_arrival = p1 * p2
                
                # --- Simulate sales for this arrival scenario (d1, d2) ---

                # 1. Class 2 customers purchase cheap tickets (priority)
                c2_sold_cheap = min(d2, c)
                c_after_c2_cheap = c - c2_sold_cheap
                d2_rem_for_exp = d2 - c2_sold_cheap

                # 2. Class 1 customers purchase cheap tickets
                c1_sold_cheap = min(d1, c_after_c2_cheap)
                c_final_after_cheap_sales = c_after_c2_cheap - c1_sold_cheap
                
                revenue_from_cheap_sales = (c2_sold_cheap + c1_sold_cheap) * PRICE_CHEAP

                # 3. Class 2 customers consider expensive tickets
                # This part is stochastic, so we calculate the expected value over its outcomes.
                expected_value_from_c2_expensive_and_future = 0
                
                # Iterate over the number of C2 customers (k) who decide to buy an expensive ticket
                for k in range(d2_rem_for_exp + 1):
                    prob_k_buys = binom_pmf(k, d2_rem_for_exp, PROB_C2_BUY_EXPENSIVE)
                    
                    # Number of C2 customers who actually purchase an expensive ticket
                    c2_sold_exp = min(e, k)
                    e_final = e - c2_sold_exp
                    
                    revenue_from_expensive_sales = c2_sold_exp * PRICE_EXPENSIVE
                    
                    # Expected revenue from the next day onwards
                    future_revenue = solve(t + 1, c_final_after_cheap_sales, e_final)
                    
                    # Accumulate the expected value for this outcome
                    expected_value_from_c2_expensive_and_future += prob_k_buys * (revenue_from_expensive_sales + future_revenue)
                
                # Total expected value for this (d1, d2) arrival scenario
                total_expected_value += prob_arrival * (revenue_from_cheap_sales + expected_value_from_c2_expensive_and_future)
            
        return total_expected_value

    # --- Find the optimal policy and maximum revenue ---
    max_revenue = -1
    best_policy_climit = -1
    results = []

    # Iterate through all possible policies (C_limit = 0 to 10)
    for c_limit in range(CAPACITY + 1):
        e_limit = CAPACITY - c_limit
        current_expected_revenue = solve(1, c_limit, e_limit)
        results.append((c_limit, current_expected_revenue))
        
        if current_expected_revenue > max_revenue:
            max_revenue = current_expected_revenue
            best_policy_climit = c_limit

    # --- Print the results ---
    print("Calculating expected revenue for different booking policies...")
    print("A policy is defined by the total number of cheap tickets offered (C_limit).")
    print("-" * 75)
    
    # This loop outputs the "numbers in the final equation"
    for c_limit, revenue in results:
        e_limit = CAPACITY - c_limit
        print(f"Policy: {c_limit:2d} cheap tickets, {e_limit:2d} expensive tickets. Expected Revenue: {revenue:.2f}")

    print("-" * 75)
    print(f"The optimal policy is to offer {best_policy_climit} cheap tickets.")
    print(f"The maximum expected total revenue is {max_revenue:.2f}.")
    
    # --- Final Answer in requested format ---
    print(f'<<<{max_revenue:.2f}>>>')

# Execute the solution
solve_revenue_maximization()