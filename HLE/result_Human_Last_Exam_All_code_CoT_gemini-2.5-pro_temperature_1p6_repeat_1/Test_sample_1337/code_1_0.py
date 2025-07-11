import collections
import sys

# To prevent deep recursion errors for the memoized recursion
sys.setrecursionlimit(2000)

def main():
    """
    Main function to find the maximum expected revenue by testing all possible booking policies (k).
    It prints the final maximum revenue found.
    """
    
    # --- Constants from the problem description ---
    TOTAL_SEATS = 10
    DAYS = 14
    C1_PRICE = 100
    C2_PRICE = 200
    C2_PURCHASE_PROB = 0.5

    # --- Customer arrival distributions ---
    N1_DIST = {0: 0.25, 1: 0.5, 2: 0.25}
    N2_DIST_WEEK1 = {0: 1.0}
    N2_DIST_WEEK2 = {0: 0.25, 1: 0.5, 2: 0.25}

    # --- Memoization table for the recursive helper ---
    # Using a nonlocal variable within the main function's scope to avoid global state.
    memo_table = {}

    def process_arrivals_recursive(s, n2, n1, c_limit):
        """
        Recursive helper to calculate expected revenue and state transitions for a specific scenario.
        It processes one customer at a time, respecting priorities and sales logic.

        Args:
            s (int): Current number of seats sold.
            n2 (int): Remaining class 2 customers to process.
            n1 (int): Remaining class 1 customers to process.
            c_limit (int): The maximum number of cheap tickets (booking limit).

        Returns:
            tuple: (expected_revenue_from_this_state, {final_s: probability_dist})
        """
        nonlocal memo_table
        state = (s, n2, n1)
        if state in memo_table:
            return memo_table[state]

        # Base case: All customers processed or all seats sold.
        if s >= TOTAL_SEATS or (n1 == 0 and n2 == 0):
            return 0.0, {s: 1.0}

        # Process Class 2 customers first due to priority for cheap tickets.
        if n2 > 0:
            if s < c_limit:
                # Sell a cheap ticket to a Class 2 customer.
                revenue, dist = process_arrivals_recursive(s + 1, n2 - 1, n1, c_limit)
                result = (C1_PRICE + revenue, dist)
            else: # Offer an expensive ticket.
                # Outcome 1: Customer buys (with C2_PURCHASE_PROB)
                rev_buy, dist_buy = process_arrivals_recursive(s + 1, n2 - 1, n1, c_limit)
                
                # Outcome 2: Customer walks away
                rev_walk, dist_walk = process_arrivals_recursive(s, n2 - 1, n1, c_limit)
                
                expected_revenue = C2_PURCHASE_PROB * (C2_PRICE + rev_buy) + (1 - C2_PURCHASE_PROB) * rev_walk
                
                final_dist = collections.defaultdict(float)
                for final_s, p in dist_buy.items():
                    final_dist[final_s] += C2_PURCHASE_PROB * p
                for final_s, p in dist_walk.items():
                    final_dist[final_s] += (1 - C2_PURCHASE_PROB) * p
                
                result = (expected_revenue, dict(final_dist))
        # Process Class 1 customers if no Class 2 customers are left.
        elif n1 > 0:
            if s < c_limit:
                # Sell a cheap ticket to a Class 1 customer.
                revenue, dist = process_arrivals_recursive(s + 1, n2, n1 - 1, c_limit)
                result = (C1_PRICE + revenue, dist)
            else:
                # Class 1 customer walks away if no cheap tickets are available.
                revenue, dist = process_arrivals_recursive(s, n2, n1 - 1, c_limit)
                result = (revenue, dist)
        else: # Should be unreachable due to the base cases
            result = (0.0, {s: 1.0})

        memo_table[state] = result
        return result

    def calculate_total_expected_revenue(k):
        """
        Calculates the total expected revenue over 14 days for a given protection level `k`.
        `k` is the number of seats protected for expensive fares.
        """
        nonlocal memo_table
        memo_table.clear()
        
        c_limit = TOTAL_SEATS - k
        
        p_state = [0.0] * (TOTAL_SEATS + 1)
        p_state[0] = 1.0
        
        total_expected_revenue = 0.0
        
        for day in range(1, DAYS + 1):
            n2_dist = N2_DIST_WEEK1 if day <= 7 else N2_DIST_WEEK2
            p_next_day_state = [0.0] * (TOTAL_SEATS + 1)
            
            for s_initial in range(TOTAL_SEATS + 1):
                if p_state[s_initial] == 0:
                    continue
                
                day_expected_revenue_from_s = 0.0
                day_final_dist_from_s = collections.defaultdict(float)

                for n1, p1 in N1_DIST.items():
                    for n2, p2 in n2_dist.items():
                        prob_scenario = p1 * p2
                        if prob_scenario == 0: continue
                        
                        rev_scenario, dist_scenario = process_arrivals_recursive(s_initial, n2, n1, c_limit)
                        
                        day_expected_revenue_from_s += prob_scenario * rev_scenario
                        for s_final, p in dist_scenario.items():
                            day_final_dist_from_s[s_final] += prob_scenario * p

                total_expected_revenue += p_state[s_initial] * day_expected_revenue_from_s
                
                for s_final, p in day_final_dist_from_s.items():
                     if s_final <= TOTAL_SEATS:
                        p_next_day_state[s_final] += p_state[s_initial] * p

            p_state = p_next_day_state
            
        return total_expected_revenue

    # --- Main execution logic ---
    max_revenue = -1.0
    
    # k is the number of seats protected for expensive fares
    for k in range(TOTAL_SEATS + 1):
        revenue_for_k = calculate_total_expected_revenue(k)
        if revenue_for_k > max_revenue:
            max_revenue = revenue_for_k

    # As per instruction, print the final answer number.
    print(f"{max_revenue:.2f}")


if __name__ == '__main__':
    main()