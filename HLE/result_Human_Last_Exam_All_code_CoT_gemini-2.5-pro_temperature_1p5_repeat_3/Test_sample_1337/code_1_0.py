import collections

# --- Problem Parameters ---
SEAT_CAPACITY = 10
TOTAL_DAYS = 14
CHEAP_PRICE = 100
EXPENSIVE_PRICE = 200

# Demand Distributions
C1_DEMAND_PROB = {0: 0.25, 1: 0.5, 2: 0.25}
C2_DEMAND_PROB = {0: 0.25, 1: 0.5, 2: 0.25}
NO_C2_DEMAND_PROB = {0: 1.0, 1: 0.0, 2: 0.0}

# --- DP State Initialization ---
# V[s] will store a tuple: (max_expected_revenue, E[cheap_sold], E[expensive_sold])
# for a state with s seats remaining.
# We iterate backwards in time, so V always holds the values for the *next* day.
# Initially, V represents values for day 15 (after all sales), which are all zero.
V = collections.defaultdict(lambda: (0.0, 0.0, 0.0))

# Memoization cache for the recursive helper function
memo = {}

def calculate_exp_value_from_exp_sales(num_c2_exp_demand, seats_available, v_next_day):
    """
    Recursively calculates the expected total value (revenue and ticket counts)
    from the point where class 2 customers are considering expensive tickets.
    This uses memoization to optimize calculations within a single day.
    """
    state = (num_c2_exp_demand, seats_available)
    if state in memo:
        return memo[state]

    # Base case: No more customers for expensive tickets today, or no seats left.
    # The only remaining value is from the future days, looked up in v_next_day.
    if num_c2_exp_demand == 0 or seats_available == 0:
        return v_next_day[seats_available]

    # --- Recursive step: consider one class 2 customer's decision ---
    
    # Case 1: Customer buys an expensive ticket (50% probability)
    rev_buy, cheap_buy, exp_buy = calculate_exp_value_from_exp_sales(
        num_c2_exp_demand - 1, seats_available - 1, v_next_day
    )
    total_rev_if_buy = EXPENSIVE_PRICE + rev_buy
    total_exp_if_buy = 1 + exp_buy
    total_cheap_if_buy = cheap_buy

    # Case 2: Customer walks away (50% probability)
    rev_walk, cheap_walk, exp_walk = calculate_exp_value_from_exp_sales(
        num_c2_exp_demand - 1, seats_available, v_next_day
    )
    
    # Calculate the expected value by averaging the two cases.
    expected_rev = 0.5 * total_rev_if_buy + 0.5 * rev_walk
    expected_cheap = 0.5 * total_cheap_if_buy + 0.5 * cheap_walk
    expected_exp = 0.5 * total_exp_if_buy + 0.5 * exp_walk

    result = (expected_rev, expected_cheap, expected_exp)
    memo[state] = result
    return result

# --- Main DP loop ---
# Iterate backwards from the last day of sales (day 14) to the first (day 1).
for day in range(TOTAL_DAYS, 0, -1):
    memo.clear()
    
    # Class 2 customers only appear in the second week (days 8-14).
    c2_demand_dist = C2_DEMAND_PROB if day > 7 else NO_C2_DEMAND_PROB
    
    V_next_day = V
    V_current_day = collections.defaultdict(lambda: (0.0, 0.0, 0.0))

    # For each number of available seats `s` at the start of the day...
    for seats in range(SEAT_CAPACITY + 1):
        best_policy_value = (-1.0, 0.0, 0.0)

        # ...find the optimal protection level `p`.
        for protection_level in range(seats + 1):
            policy_exp_rev, policy_exp_cheap, policy_exp_exp = 0.0, 0.0, 0.0

            # Iterate over all possible demand scenarios for C1 and C2.
            for c1_demand, c1_prob in C1_DEMAND_PROB.items():
                for c2_demand, c2_prob in c2_demand_dist.items():
                    prob = c1_prob * c2_prob
                    if prob == 0:
                        continue

                    # Number of cheap tickets we are willing to sell today.
                    cheap_tickets_limit = seats - protection_level

                    # C2 customers have priority for cheap tickets.
                    c2_bought_cheap = min(c2_demand, cheap_tickets_limit)

                    # C1 customers buy remaining cheap tickets.
                    c1_bought_cheap = min(c1_demand, cheap_tickets_limit - c2_bought_cheap)

                    cheap_sold_today = c1_bought_cheap + c2_bought_cheap
                    rev_from_cheap_today = cheap_sold_today * CHEAP_PRICE

                    seats_after_cheap = seats - cheap_sold_today
                    c2_demand_for_expensive = c2_demand - c2_bought_cheap

                    # Calculate expected value from expensive sales today + all future sales.
                    (exp_future_rev, 
                     exp_future_cheap, 
                     exp_future_exp) = calculate_exp_value_from_exp_sales(
                        c2_demand_for_expensive, seats_after_cheap, V_next_day
                    )
                    
                    scenario_total_rev = rev_from_cheap_today + exp_future_rev
                    scenario_total_cheap = cheap_sold_today + exp_future_cheap
                    scenario_total_exp = exp_future_exp

                    policy_exp_rev += prob * scenario_total_rev
                    policy_exp_cheap += prob * scenario_total_cheap
                    policy_exp_exp += prob * scenario_total_exp
            
            # If this policy is better, update the best value found so far.
            if policy_exp_rev > best_policy_value[0]:
                best_policy_value = (policy_exp_rev, policy_exp_cheap, policy_exp_exp)

        # Store the optimal value for the state (day, seats).
        V_current_day[seats] = best_policy_value

    V = V_current_day

# The final result is the optimal value for Day 1 with full capacity.
final_revenue, final_exp_cheap, final_exp_exp = V[SEAT_CAPACITY]

print("Final Equation for Maximum Expected Revenue:")
print("E[Cheap Tickets Sold] * Price + E[Expensive Tickets Sold] * Price = E[Total Revenue]")
print(f"{final_exp_cheap:.4f} * {CHEAP_PRICE} + {final_exp_exp:.4f} * {EXPENSIVE_PRICE} = {final_revenue:.4f}")
