import math

# Use a global dictionary for memoization to store results of subproblems
memo = {}

# --- Problem Parameters ---
CAPACITY = 10
TOTAL_DAYS = 14
CHEAP_PRICE = 100
EXPENSIVE_PRICE = 200

# Customer arrival distributions
N1_DIST = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1 customers (all days)
N2_DIST_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2 customers (days 8-14)
N2_DIST_WEEK1 = {0: 1.0}  # Class 2 customers (days 1-7 have no arrivals)


def get_expected_revenue(day, seats_rem, cheap_sold, q):
    """
    Calculates the maximum expected future revenue using dynamic programming.
    The state is defined by (day, seats_rem, cheap_sold).
    'q' is the fixed protection level for this entire calculation.
    """
    # Base case: If selling period is over or no seats are left, no more revenue.
    if day > TOTAL_DAYS or seats_rem == 0:
        return 0

    # Memoization: Return stored result if subproblem was already solved.
    state = (day, seats_rem, cheap_sold)
    if state in memo:
        return memo[state]

    # Determine the arrival distribution for Class 2 customers for the current day.
    n2_dist = N2_DIST_WEEK1 if day <= 7 else N2_DIST_WEEK2

    total_expected_value = 0
    # Iterate over all possible numbers of arriving customers (n1, n2)
    # weighted by their probabilities (p1, p2).
    for n1, p1 in N1_DIST.items():
        for n2, p2 in n2_dist.items():

            exp_val_for_n1n2 = 0
            cheap_ticket_limit = CAPACITY - q

            # Step 1: Process Class 2 Customers (they have priority)
            cheap_available = max(0, cheap_ticket_limit - cheap_sold)
            c2_buy_cheap = min(n2, cheap_available, seats_rem)
            rev_c2_cheap = c2_buy_cheap * CHEAP_PRICE
            
            s_after_c2_cheap = seats_rem - c2_buy_cheap
            c_after_c2_cheap = cheap_sold + c2_buy_cheap

            # Remaining Class 2 customers are offered expensive tickets.
            c2_offered_exp = min(n2 - c2_buy_cheap, s_after_c2_cheap)
            
            # The number of expensive sales 'b' is probabilistic.
            # We sum over all possible outcomes for 'b' using the binomial distribution.
            # This loop correctly handles the c2_offered_exp = 0 case.
            prob_factor_exp = 0.5 ** c2_offered_exp
            for b in range(c2_offered_exp + 1):
                prob_b = math.comb(c2_offered_exp, b) * prob_factor_exp
                
                rev_from_exp_sales = b * EXPENSIVE_PRICE
                s_after_exp_sales = s_after_c2_cheap - b
                c_after_exp_sales = c_after_c2_cheap
                
                # Step 2: Process Class 1 Customers
                cheap_avail_for_c1 = max(0, cheap_ticket_limit - c_after_exp_sales)
                c1_buy_cheap = min(n1, cheap_avail_for_c1, s_after_exp_sales)
                rev_from_c1_sales = c1_buy_cheap * CHEAP_PRICE

                day_revenue = rev_c2_cheap + rev_from_exp_sales + rev_from_c1_sales
                
                # Determine the state for the next day.
                s_next = s_after_exp_sales - c1_buy_cheap
                c_next = c_after_exp_sales + c1_buy_cheap
                
                # Recursively call for the expected future revenue.
                future_revenue = get_expected_revenue(day + 1, s_next, c_next, q)
                
                exp_val_for_n1n2 += prob_b * (day_revenue + future_revenue)

            total_expected_value += p1 * p2 * exp_val_for_n1n2
            
    # Store the result in the memoization table and return it.
    memo[state] = total_expected_value
    return total_expected_value

def find_max_revenue():
    """
    Iterates through all possible protection levels (q) to find the optimal
    policy and the corresponding maximum expected revenue.
    """
    max_revenue = 0.0
    optimal_q = -1

    # Loop through all possible protection levels q (from 0 to capacity).
    for q in range(CAPACITY + 1):
        global memo
        memo = {}
        # Calculate expected revenue for the entire period with policy 'q'.
        revenue_for_q = get_expected_revenue(1, CAPACITY, 0, q)
        
        if revenue_for_q > max_revenue:
            max_revenue = revenue_for_q
            optimal_q = q
            
    print("To find the maximum expected revenue, we test every possible booking policy.")
    print("A policy is defined by a 'protection level' (q), which is the number of seats (out of 10) reserved for high-fare tickets.")
    print(f"\nThe optimal policy is to have a protection level of q = {optimal_q}.")
    print(f"This means reserving {optimal_q} seats for expensive tickets (${EXPENSIVE_PRICE}) and allowing a maximum of {CAPACITY - optimal_q} cheap tickets (${CHEAP_PRICE}) to be sold.")
    print("\nThe final calculation for the maximum expected revenue is:")
    print(f"Max E[Revenue] = E[Revenue(q={optimal_q})] = ${max_revenue:.2f}")

if __name__ == '__main__':
    find_max_revenue()