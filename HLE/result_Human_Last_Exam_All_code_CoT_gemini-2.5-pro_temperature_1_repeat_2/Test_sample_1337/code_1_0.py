import collections
from math import comb

# --- Helper function for binomial probability ---
# Using a cache to speed up repeated calculations of the same binomial probability.
memo_binomial = {}
def binomial_pmf(k, n, p):
    """
    Calculates the probability mass function of a binomial distribution B(n, p) at k.
    """
    if (k, n, p) in memo_binomial:
        return memo_binomial[(k, n, p)]
    if k < 0 or k > n:
        return 0
    # Using the combination formula: C(n, k) * p^k * (1-p)^(n-k)
    result = comb(n, k) * (p**k) * ((1-p)**(n-k))
    memo_binomial[(k, n, p)] = result
    return result

def solve_revenue_maximization():
    """
    Solves the airline revenue maximization problem using dynamic programming.
    """
    # --- Problem Parameters ---
    SEAT_CAPACITY = 10
    DAYS = 14
    CHEAP_PRICE = 100
    EXPENSIVE_PRICE = 200

    # Probability distribution for Class 1 customer arrivals (same for all 14 days)
    PROB_D1 = {0: 0.25, 1: 0.5, 2: 0.25}
    # Probability distribution for Class 2 customer arrivals
    PROB_D2_WEEK1 = {0: 1.0}  # Days 1-7: No Class 2 customers
    PROB_D2_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25} # Days 8-14

    # --- Dynamic Programming Setup ---
    # DP table: V[t][s] stores the maximum expected revenue from day t onwards
    # with s seats remaining.
    # Dimensions: (DAYS + 2) rows for days 1 to 14 (and base cases t=15, t=0),
    #             (SEAT_CAPACITY + 1) columns for seats 0 to 10.
    V = [[0.0 for _ in range(SEAT_CAPACITY + 1)] for _ in range(DAYS + 2)]
    
    # We iterate backwards in time, from the last day of sales (t=14) to the first (t=1).
    # The base case V[15][s] = 0 is already initialized (no more revenue after day 14).
    for t in range(DAYS, 0, -1):
        
        # Determine the Class 2 demand distribution for the current day t.
        # Days 1-7 are the first week of sales, Days 8-14 are the second week.
        if 1 <= t <= 7:
            prob_d2_dist = PROB_D2_WEEK1
        else:  # 8 <= t <= 14
            prob_d2_dist = PROB_D2_WEEK2
        
        # Iterate over all possible numbers of remaining seats `s` at the start of day `t`.
        for s in range(SEAT_CAPACITY + 1):
            max_expected_revenue_for_s = -1.0
            
            # The decision variable is `k`: the number of cheap tickets to make available.
            # We test every possible value of `k` from 0 to `s`.
            for k in range(s + 1):
                current_k_expected_revenue = 0.0
                
                # To calculate the expected revenue for a given `k`, we sum over all
                # possible demand scenarios (d1, d2).
                for d1, p1 in PROB_D1.items():
                    for d2, p2 in prob_d2_dist.items():
                        prob_demand = p1 * p2
                        
                        # --- Sales Process Simulation for a given (d1, d2) demand ---

                        # 1. Class 2 customers have priority for cheap tickets.
                        sold_c2_cheap = min(d2, k)
                        revenue_from_cheap_tickets = CHEAP_PRICE * sold_c2_cheap
                        
                        remaining_k_for_c1 = k - sold_c2_cheap
                        remaining_d2_for_exp = d2 - sold_c2_cheap

                        # 2. Class 1 customers buy the remaining cheap tickets.
                        sold_c1_cheap = min(d1, remaining_k_for_c1)
                        revenue_from_cheap_tickets += CHEAP_PRICE * sold_c1_cheap
                        
                        # 3. Remaining Class 2 customers consider expensive tickets.
                        # This part involves an expectation over the 50% chance of buying.
                        num_expensive_seats_available = s - k
                        exp_value_from_exp_and_future = 0.0

                        if remaining_d2_for_exp > 0 and num_expensive_seats_available > 0:
                            # We calculate the expected value by summing over all possible outcomes
                            # of how many C2 customers decide to buy an expensive ticket.
                            for i in range(remaining_d2_for_exp + 1):
                                prob_i_buy = binomial_pmf(i, remaining_d2_for_exp, 0.5)
                                
                                sold_c2_exp = min(i, num_expensive_seats_available)
                                revenue_from_expensive_tickets = EXPENSIVE_PRICE * sold_c2_exp
                                
                                seats_sold_total = sold_c1_cheap + sold_c2_cheap + sold_c2_exp
                                seats_left = s - seats_sold_total
                                
                                # Add the future revenue from the next state V[t+1][seats_left]
                                future_revenue = V[t + 1][seats_left]
                                
                                exp_value_from_exp_and_future += prob_i_buy * (revenue_from_expensive_tickets + future_revenue)
                        else:
                            # If no C2 customers are left for expensive tickets or none are available,
                            # the only future revenue is from the state determined by cheap sales.
                            seats_sold_total = sold_c1_cheap + sold_c2_cheap
                            seats_left = s - seats_sold_total
                            future_revenue = V[t + 1][seats_left]
                            exp_value_from_exp_and_future = future_revenue
                        
                        # Add the weighted value of this scenario to the total expected revenue for this `k`.
                        current_k_expected_revenue += prob_demand * (revenue_from_cheap_tickets + exp_value_from_exp_and_future)
                
                # Update the maximum revenue for state (t, s) if this `k` is better.
                if current_k_expected_revenue > max_expected_revenue_for_s:
                    max_expected_revenue_for_s = current_k_expected_revenue
            
            V[t][s] = max_expected_revenue_for_s
            
    # The final answer is the maximum expected revenue at the beginning,
    # on Day 1 with the full seat capacity.
    max_revenue = V[1][SEAT_CAPACITY]
    
    # --- Output the final result as an equation ---
    print(f"To find the maximum expected revenue, a dynamic programming model is used.")
    print(f"The calculation finds the optimal booking policy for each day based on the number of remaining seats.")
    print(f"Based on the parameters:\n- Seat Capacity: {SEAT_CAPACITY}\n- Days: {DAYS}\n- Cheap Price: {CHEAP_PRICE}\n- Expensive Price: {EXPENSIVE_PRICE}")
    print("\nThe final equation for the maximum expected total revenue is:")
    # The following line presents the final answer in an equation-like format
    # summarizing the inputs and the final calculated value.
    print(f"MaxRevenue(seats={SEAT_CAPACITY}, days={DAYS}, prices=({CHEAP_PRICE}, {EXPENSIVE_PRICE})) = {max_revenue:.2f}")
    print(f"<<<{max_revenue:.2f}>>>")

# Execute the solver
solve_revenue_maximization()