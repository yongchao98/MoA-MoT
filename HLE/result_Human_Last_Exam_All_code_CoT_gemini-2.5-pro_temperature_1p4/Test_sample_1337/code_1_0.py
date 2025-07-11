import math

def solve_revenue_maximization():
    """
    Solves the airline revenue management problem using dynamic programming
    to find the maximum expected total revenue.
    """
    # --- Problem Constants ---
    TOTAL_SEATS = 10
    TOTAL_DAYS = 14
    PRICE_CHEAP = 100
    PRICE_EXPENSIVE = 200

    # --- Customer Arrival Distributions ---
    # {number_of_requests: probability}
    CLASS_1_DIST = {0: 0.25, 1: 0.50, 2: 0.25}
    CLASS_2_DIST = {0: 0.25, 1: 0.50, 2: 0.25}

    def get_prob_k_sales(k, num_customers, num_seats):
        """
        Calculates the probability of selling exactly k expensive tickets.
        Args:
            k (int): The number of sales.
            num_customers (int): The number of customers considering an expensive ticket.
            num_seats (int): The number of available seats.
        """
        if num_customers == 0:
            return 1.0 if k == 0 else 0.0
        if num_seats == 0:
            return 1.0 if k == 0 else 0.0

        if k > min(num_customers, num_seats):
            return 0.0
        
        # Pre-calculated probabilities for N=1 and N=2
        if num_customers == 1:
            if k == 0: return 0.5
            if k == 1: return 0.5
        elif num_customers == 2:
            if num_seats == 1:
                if k == 0: return 0.25
                if k == 1: return 0.75
            else:  # num_seats >= 2
                if k == 0: return 0.25
                if k == 1: return 0.50
                if k == 2: return 0.25
        return 0.0

    # --- Dynamic Programming Setup ---
    # V[t][s]: max expected revenue from day t to 14 with s seats left.
    # We use t from 1 to 15. V[15] is the base case (revenue=0 after day 14).
    V = [[0.0] * (TOTAL_SEATS + 1) for _ in range(TOTAL_DAYS + 2)]

    # --- DP Calculation (Iterate backwards in time) ---
    for t in range(TOTAL_DAYS, 0, -1):
        # Determine Class 2 customer distribution for the current day
        current_c2_dist = CLASS_2_DIST if t > 7 else {0: 1.0}

        # Iterate over all possible numbers of remaining seats
        for s in range(TOTAL_SEATS + 1):
            if s == 0:
                V[t][s] = 0
                continue
            
            max_expected_revenue_for_s = -1.0

            # Iterate over all possible decisions (protection level b)
            for b in range(s + 1):
                cheap_ticket_limit = s - b
                current_b_expected_revenue = 0.0
                
                # Sum over all possible customer arrival scenarios (n1, n2)
                for n1, p1 in CLASS_1_DIST.items():
                    for n2, p2 in current_c2_dist.items():
                        prob_scenario = p1 * p2
                        
                        # 1. Process cheap ticket sales (C2 has priority)
                        c2_buy_cheap = min(n2, cheap_ticket_limit)
                        c1_buy_cheap = min(n1, cheap_ticket_limit - c2_buy_cheap)
                        
                        revenue_from_cheap = (c2_buy_cheap + c1_buy_cheap) * PRICE_CHEAP
                        seats_after_cheap = s - (c2_buy_cheap + c1_buy_cheap)
                        n2_rem = n2 - c2_buy_cheap

                        # 2. Process expensive ticket sales (probabilistic)
                        # We calculate the expected value from this stage onwards.
                        expected_value_from_expensive_stage = 0.0
                        
                        # Iterate over possible numbers of expensive sales (k)
                        for k in range(min(n2_rem, seats_after_cheap) + 1):
                            prob_k = get_prob_k_sales(k, n2_rem, seats_after_cheap)
                            
                            revenue_from_expensive = k * PRICE_EXPENSIVE
                            final_seats = seats_after_cheap - k
                            
                            # V[t+1] holds the max expected future revenue from the next day
                            future_revenue = V[t + 1][final_seats]
                            
                            expected_value_from_expensive_stage += prob_k * (revenue_from_expensive + future_revenue)
                            
                        current_b_expected_revenue += prob_scenario * (revenue_from_cheap + expected_value_from_expensive_stage)

                if current_b_expected_revenue > max_expected_revenue_for_s:
                    max_expected_revenue_for_s = current_b_expected_revenue
            
            V[t][s] = max_expected_revenue_for_s

    # The final answer is the max expected revenue at the start (Day 1, 10 seats)
    final_revenue = V[1][TOTAL_SEATS]
    
    # As requested, here is the equation for the final answer.
    # The final answer is V(day=1, seats=10), which is the result of the entire DP calculation.
    print(f"Maximum Expected Total Revenue = {final_revenue:.4f}")

solve_revenue_maximization()