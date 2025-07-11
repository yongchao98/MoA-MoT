import math
from functools import lru_cache

# Using lru_cache for automatic memoization, which is cleaner than a manual dictionary.
# maxsize is set to None for an unbounded cache.
@lru_cache(maxsize=None)
def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

@lru_cache(maxsize=None)
def expected_values(d, c_rem, e_rem):
    """
    Calculates expected future values using dynamic programming.
    
    Args:
        d (int): The current day (from 1 to 14).
        c_rem (int): Number of cheap tickets remaining.
        e_rem (int): Number of expensive tickets remaining.
        
    Returns:
        tuple: (expected_revenue, expected_cheap_sold, expected_expensive_sold)
    """
    # Base case: End of selling period or no tickets left.
    if d > 14 or (c_rem <= 0 and e_rem <= 0):
        return (0, 0, 0)

    # Define daily demand distributions
    n1_dist = {0: 0.25, 1: 0.5, 2: 0.25}
    n2_dist = {0: 1.0} if d <= 7 else {0: 0.25, 1: 0.5, 2: 0.25}

    total_exp_rev = 0
    total_exp_cheap_sold = 0
    total_exp_exp_sold = 0

    # Iterate over all possible demand scenarios for the day
    for n1, p1 in n1_dist.items():
        for n2, p2 in n2_dist.items():
            prob_scenario = p1 * p2
            if prob_scenario == 0:
                continue

            # --- Sales process for this (n1, n2) scenario ---
            
            # 1. Class 2 customers buy cheap tickets first (priority)
            c_sold_to_2 = min(n2, c_rem)
            c_rem_after_c2 = c_rem - c_sold_to_2
            
            # 2. Class 1 customers buy remaining cheap tickets
            c_sold_to_1 = min(n1, c_rem_after_c2)
            c_rem_final = c_rem_after_c2 - c_sold_to_1
            
            c_sold_today = c_sold_to_2 + c_sold_to_1

            # 3. Unhappy Class 2 customers consider expensive tickets
            n2_unhappy = n2 - c_sold_to_2
            
            # These are the expected values *given* the (n1, n2) scenario
            exp_rev_scenario = 0
            exp_cheap_sold_scenario = 0
            exp_exp_sold_scenario = 0

            if n2_unhappy > 0 and e_rem > 0:
                # Average over the binomial outcomes for customers wanting an expensive ticket
                for k in range(n2_unhappy + 1):
                    prob_k = combinations(n2_unhappy, k) * (0.5)**n2_unhappy
                    
                    e_sold_k = min(k, e_rem)
                    e_rem_final_k = e_rem - e_sold_k
                    
                    future_rev, future_cheap, future_exp = expected_values(d + 1, c_rem_final, e_rem_final_k)
                    
                    # Add the contribution of this outcome k to the scenario's expectation
                    exp_rev_scenario += prob_k * ((c_sold_today * 100 + e_sold_k * 200) + future_rev)
                    exp_cheap_sold_scenario += prob_k * (c_sold_today + future_cheap)
                    exp_exp_sold_scenario += prob_k * (e_sold_k + future_exp)
            else:
                # Deterministic case (no unhappy C2 or no expensive tickets)
                future_rev, future_cheap, future_exp = expected_values(d + 1, c_rem_final, e_rem)
                exp_rev_scenario = (c_sold_today * 100) + future_rev
                exp_cheap_sold_scenario = c_sold_today + future_cheap
                exp_exp_sold_scenario = 0 + future_exp

            # Add the scenario's expected values to the daily total, weighted by probability
            total_exp_rev += prob_scenario * exp_rev_scenario
            total_exp_cheap_sold += prob_scenario * exp_cheap_sold_scenario
            total_exp_exp_sold += prob_scenario * exp_exp_sold_scenario

    return (total_exp_rev, total_exp_cheap_sold, total_exp_exp_sold)

def find_max_revenue():
    """
    Iterates through all policies to find the one with the maximum expected revenue.
    """
    best_policy = {'p': -1, 'revenue': -1, 'cheap_sold': 0, 'exp_sold': 0}

    # A policy is defined by the protection level 'p'
    for p in range(11):  # p from 0 to 10
        c_initial = 10 - p
        e_initial = p
        
        # Clear cache for each policy to ensure independent calculation
        expected_values.cache_clear()
        
        rev, cheap_sold, exp_sold = expected_values(1, c_initial, e_initial)
        
        if rev > best_policy['revenue']:
            best_policy['p'] = p
            best_policy['revenue'] = rev
            best_policy['cheap_sold'] = cheap_sold
            best_policy['exp_sold'] = exp_sold
            
    # Output the results for the best policy found
    p = best_policy['p']
    num_cheap = 10 - p
    num_expensive = p
    max_revenue = best_policy['revenue']
    e_cheap_sold = best_policy['cheap_sold']
    e_exp_sold = best_policy['exp_sold']

    print(f"The optimal policy is to set a booking limit of {num_cheap} for cheap tickets and protect {num_expensive} seats for expensive fares.")
    print("\nThe final equation for the maximum expected revenue is:")
    print(f"E[Revenue] = E[Cheap Sold] * $100 + E[Expensive Sold] * $200")
    print(f"{max_revenue:.2f} = {e_cheap_sold:.4f} * 100 + {e_exp_sold:.4f} * 200")

    # Final answer in the required format
    print(f"\nMaximum Expected Total Revenue:")
    print(f"{max_revenue:.2f}")
    
    return max_revenue

if __name__ == '__main__':
    max_rev = find_max_revenue()
    print(f"<<<{max_rev:.2f}>>>")