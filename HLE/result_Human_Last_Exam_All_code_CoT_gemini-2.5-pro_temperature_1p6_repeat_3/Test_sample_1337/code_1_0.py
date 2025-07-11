import math
from functools import lru_cache

# --- Problem Parameters ---
C1_DEMAND_DIST = {0: 0.25, 1: 0.50, 2: 0.25}
C2_DEMAND_DIST = {0: 0.25, 1: 0.50, 2: 0.25}
C2_NO_DEMAND_DIST = {0: 1.0}
CHEAP_PRICE = 100
EXPENSIVE_PRICE = 200
TOTAL_SEATS = 10
TOTAL_DAYS = 14
PROB_C2_BUY_EXPENSIVE = 0.5

@lru_cache(maxsize=None)
def get_sales_dist(num_requests, num_items, prob_buy):
    """
    Calculates the probability distribution of the number of items sold given
    probabilistic buyers.
    Args:
        num_requests: The number of customers.
        num_items: The number of items available.
        prob_buy: The probability that a single customer will make a purchase.
    Returns:
        A dictionary {num_sales: probability}.
    """
    if num_requests == 0 or num_items == 0:
        return {0: 1.0}

    # prob_k_wants[k] is the probability that exactly k customers want to buy
    prob_k_wants = {
        k: math.comb(num_requests, k) * (prob_buy**k) * ((1 - prob_buy)**(num_requests - k))
        for k in range(num_requests + 1)
    }

    dist = {}
    # For sales < num_items, it means exactly that many customers wanted to buy
    for k_sales in range(num_items):
        dist[k_sales] = prob_k_wants.get(k_sales, 0)
    
    # For sales == num_items, it means num_items or more customers wanted to buy
    prob_at_least_k_items = sum(
        prob_k_wants.get(k, 0) for k in range(num_items, num_requests + 1)
    )
    dist[num_items] = prob_at_least_k_items
    
    return dist

@lru_cache(maxsize=None)
def solve_dp(day, cheap_seats_left, expensive_seats_left):
    """
    Calculates the maximum expected revenue using dynamic programming.
    The state is defined by (day, cheap_seats_left, expensive_seats_left).
    """
    if day > TOTAL_DAYS or (cheap_seats_left == 0 and expensive_seats_left == 0):
        return 0

    c2_dist = C2_DEMAND_DIST if day > 7 else C2_NO_DEMAND_DIST
    c1_dist = C1_DEMAND_DIST
    total_expected_revenue = 0

    for d1, p1 in c1_dist.items():
        for d2, p2 in c2_dist.items():
            scenario_prob = p1 * p2
            if scenario_prob == 0:
                continue

            # 1. Class 2 customers purchase cheap tickets (priority)
            c2_cheap_sales = min(d2, cheap_seats_left)
            rev_from_c2_cheap = c2_cheap_sales * CHEAP_PRICE
            remaining_cheap_after_c2 = cheap_seats_left - c2_cheap_sales
            unfulfilled_c2 = d2 - c2_cheap_sales

            # 2. Class 1 customers purchase cheap tickets
            c1_cheap_sales = min(d1, remaining_cheap_after_c2)
            rev_from_c1_cheap = c1_cheap_sales * CHEAP_PRICE
            final_cheap_seats = remaining_cheap_after_c2 - c1_cheap_sales

            # 3. Unfulfilled Class 2 customers consider expensive tickets
            expected_value_from_expensive_leg = 0
            if unfulfilled_c2 > 0 and expensive_seats_left > 0:
                sales_dist = get_sales_dist(unfulfilled_c2, expensive_seats_left, PROB_C2_BUY_EXPENSIVE)
                
                for num_sales, prob in sales_dist.items():
                    if prob == 0:
                        continue
                    rev_from_expensive = num_sales * EXPENSIVE_PRICE
                    remaining_expensive = expensive_seats_left - num_sales
                    future_revenue = solve_dp(day + 1, final_cheap_seats, remaining_expensive)
                    expected_value_from_expensive_leg += prob * (rev_from_expensive + future_revenue)
            else:
                future_revenue = solve_dp(day + 1, final_cheap_seats, expensive_seats_left)
                expected_value_from_expensive_leg = future_revenue

            scenario_total_value = rev_from_c2_cheap + rev_from_c1_cheap + expected_value_from_expensive_leg
            total_expected_revenue += scenario_prob * scenario_total_value
            
    return total_expected_revenue

def find_max_revenue():
    """
    Iterates through all possible booking policies to find the maximum expected revenue.
    """
    max_revenue = 0
    optimal_policy = (0, 0)
    revenues = []
    
    print("Calculating expected revenue for each booking policy (c, e)...")
    print("-" * 40)
    print("Policy(Cheap, Exp) | Exp. Revenue")
    print("-" * 40)

    for c in range(TOTAL_SEATS + 1):
        e = TOTAL_SEATS - c
        expected_revenue = solve_dp(1, c, e)
        revenues.append(expected_revenue)
        
        print(f"({c:2d}, {e:2d})              | ${expected_revenue:8.2f}")

        if expected_revenue > max_revenue:
            max_revenue = expected_revenue
            optimal_policy = (c, e)

    print("-" * 40)
    print(f"\nOptimal Policy: {optimal_policy[0]} cheap tickets, {optimal_policy[1]} expensive tickets.")
    
    # Per instruction, showing the numbers in the final equation
    print("\nThe final result is the maximum of the expected revenues calculated for each policy:")
    revenue_strings = [f"{r:.2f}" for r in revenues]
    print(f"Max Revenue = max([{', '.join(revenue_strings)}])")
    
    print(f"\nThe maximum expected total revenue is: ${max_revenue:.2f}")
    return max_revenue

if __name__ == "__main__":
    final_answer = find_max_revenue()
    # The final answer in the requested format
    # print(f"\n<<<{final_answer:.2f}>>>")