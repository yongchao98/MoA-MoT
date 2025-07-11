import math
import sys

# It's good practice to increase the recursion limit for deep DP problems.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems to avoid re-computation.
memo = {}

# Define the probability distributions for customer arrivals.
P_C1 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1 customers (all days)
P_C2_W1 = {0: 1.0, 1: 0.0, 2: 0.0}     # Class 2 customers (Days 1-7, Week 1)
P_C2_W2 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2 customers (Days 8-14, Week 2)

def binom_prob(n, k, p):
    """
    Calculates the binomial probability P(X=k) for a random variable X ~ B(n, p).
    This is used for the Class 2 customers' decision to buy an expensive ticket.
    """
    if k < 0 or k > n:
        return 0
    if p == 0:
        return 1.0 if k == 0 else 0.0
    if p == 1:
        return 1.0 if k == n else 0.0
    # Use math.comb for combinations (nCk)
    return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

def calculate_expected_revenue(d, n_c, n_e):
    """
    Calculates the expected future revenue using dynamic programming with memoization.
    The state is defined by (d, n_c, n_e):
      d: current day (from 1 to 14)
      n_c: number of cheap tickets ($100) remaining
      n_e: number of expensive tickets ($200) remaining
    """
    # Base case: After the selling period (day 14), no more revenue can be generated.
    if d > 14:
        return 0

    state = (d, n_c, n_e)
    if state in memo:
        return memo[state]

    # Determine Class 2 customer arrival probabilities based on the current week.
    P_C2 = P_C2_W1 if d <= 7 else P_C2_W2
    
    total_expected_revenue_for_state = 0
    
    # Iterate over all possible numbers of arriving customers for the day.
    for c1_req, p_c1 in P_C1.items():
        for c2_req, p_c2 in P_C2.items():
            
            prob_of_this_arrival_event = p_c1 * p_c2
            if prob_of_this_arrival_event == 0:
                continue

            # --- Simulate the sales process for this specific arrival event ---
            
            # 1. Class 2 customers have priority for cheap tickets.
            c2_sold_cheap = min(c2_req, n_c)
            rev_from_c2_cheap = c2_sold_cheap * 100
            
            rem_n_c_after_c2 = n_c - c2_sold_cheap
            rem_c2_customers = c2_req - c2_sold_cheap
            
            # 2. Class 1 customers buy from the remaining pool of cheap tickets.
            c1_sold_cheap = min(c1_req, rem_n_c_after_c2)
            rev_from_c1_cheap = c1_sold_cheap * 100
            
            final_n_c = rem_n_c_after_c2 - c1_sold_cheap
            day_rev_from_cheap_sales = rev_from_c1_cheap + rev_from_c2_cheap

            # 3. Remaining Class 2 customers consider buying expensive tickets.
            exp_val_from_exp_sales_and_future = 0
            
            if rem_c2_customers > 0 and n_e > 0:
                # Sum over all possibilities for how many Class 2 customers purchase an expensive ticket.
                for j in range(rem_c2_customers + 1):
                    prob_j_wants_to_buy = binom_prob(rem_c2_customers, j, 0.5)
                    
                    c2_sold_exp = min(j, n_e)
                    rev_from_c2_exp = c2_sold_exp * 200
                    
                    final_n_e = n_e - c2_sold_exp
                    
                    future_revenue = calculate_expected_revenue(d + 1, final_n_c, final_n_e)
                    exp_val_from_exp_sales_and_future += prob_j_wants_to_buy * (rev_from_c2_exp + future_revenue)
            else:
                # If no Class 2 customers can buy expensive tickets, value comes from future days.
                future_revenue = calculate_expected_revenue(d + 1, final_n_c, n_e)
                exp_val_from_exp_sales_and_future = future_revenue
                
            total_expected_revenue_for_state += prob_of_this_arrival_event * (day_rev_from_cheap_sales + exp_val_from_exp_sales_and_future)
            
    # Memoize and return the result for the current state.
    memo[state] = total_expected_revenue_for_state
    return total_expected_revenue_for_state

def find_optimal_policy_and_revenue():
    """
    Iterates through all possible booking policies to find the one that maximizes total expected revenue.
    """
    max_revenue = 0.0
    best_policy_k = -1

    # A policy is defined by 'k', the number of cheap tickets to offer (out of 10 total seats).
    for k in range(11):  # k from 0 to 10
        num_cheap_initial = k
        num_expensive_initial = 10 - k
        
        global memo
        memo.clear()
        
        # Calculate the total expected revenue for this policy, starting from day 1.
        revenue_for_policy = calculate_expected_revenue(1, num_cheap_initial, num_expensive_initial)
        
        if revenue_for_policy > max_revenue:
            max_revenue = revenue_for_policy
            best_policy_k = k

    optimal_num_cheap = best_policy_k
    optimal_num_expensive = 10 - best_policy_k
    
    print(f"Maximum Expected Revenue = Revenue(cheap_tickets={optimal_num_cheap}, expensive_tickets={optimal_num_expensive}) = {max_revenue}")
    print(f"<<<{max_revenue}>>>")

# Run the main function to find and print the result.
find_optimal_policy_and_revenue()