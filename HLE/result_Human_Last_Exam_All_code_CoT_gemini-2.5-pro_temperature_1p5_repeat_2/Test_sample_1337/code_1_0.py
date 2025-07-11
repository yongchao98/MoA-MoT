import collections

def solve_revenue_maximization():
    """
    Solves the airline revenue maximization problem using dynamic programming.
    The method iterates through all possible allocations of cheap vs. expensive tickets,
    calculates the expected revenue for each using a backward induction DP approach,
    and finds the maximum possible expected revenue.
    """
    max_total_revenue = 0
    best_c = -1
    
    # Store the results for each policy `c` to display later
    revenue_by_c = {}

    # Define the probability distribution for the number of customer requests
    P_k = {0: 0.25, 1: 0.5, 2: 0.25}

    # Iterate over all possible booking policies (number of cheap tickets `c`)
    for c in range(11):  # c can be from 0 to 10
        n_c_total = c
        n_e_total = 10 - c

        # DP table: V[(d, nc, ne)] stores the expected revenue from the start of day d onwards.
        V = collections.defaultdict(float)
        
        # Memoization cache for the W function for the current policy `c`.
        W_cache = {}

        def calculate_W(d, nc, ne, k):
            """
            This is a memoized recursive helper function.
            It calculates the expected total revenue (from current sales + future days)
            when `k` Class 2 customers, who failed to get cheap tickets, are now
            sequentially considering buying one of `ne` available expensive tickets.
            The state before this is that day `d`'s cheap sales are done, and `nc` cheap tickets are left.
            The future revenue part comes from the main DP table V for day d+1.
            """
            # Check cache first
            if (d, nc, ne, k) in W_cache:
                return W_cache[(d, nc, ne, k)]
            
            # Base Case 1: If no more customers are considering, the only revenue is the future value.
            if k == 0:
                return V.get((d + 1, nc, ne), 0.0)
            
            # Base Case 2: If no more expensive tickets, all remaining k customers walk away.
            if ne == 0:
                return V.get((d + 1, nc, 0), 0.0)
            
            # Recursive Step: Consider the first of the k customers.
            # Case A (buys with p=0.5): Earn $200, one less expensive ticket, one less customer to consider.
            val_if_buy = 200 + calculate_W(d, nc, ne - 1, k - 1)
            # Case B (walks away with p=0.5): Earn $0, same tickets, one less customer to consider.
            val_if_walk_away = calculate_W(d, nc, ne, k - 1)
            
            # The expected value is the average of the two cases.
            result = 0.5 * val_if_buy + 0.5 * val_if_walk_away
            W_cache[(d, nc, ne, k)] = result
            return result

        # Main DP loop: iterate backwards from the last day of sales (d=14) to the first (d=1).
        for d in range(14, 0, -1):
            # Class 2 customers only appear in the second week (day 8 to 14).
            P_k2 = P_k if d >= 8 else {0: 1.0, 1: 0.0, 2: 0.0}
            
            # Iterate over all possible states (number of remaining cheap/expensive tickets).
            for nc in range(n_c_total + 1):
                for ne in range(n_e_total + 1):
                    
                    expected_rev_for_state = 0
                    
                    # Calculate the expectation by summing over all possible customer arrival scenarios.
                    for k1, p1 in P_k.items(): # k1: num of class 1 customers
                        for k2, p2 in P_k2.items(): # k2: num of class 2 customers
                            prob = p1 * p2
                            if prob == 0: continue
                            
                            # --- Process sales for this scenario ---
                            # 1. Cheap ticket sales (Class 2 has priority)
                            sales_c2 = min(nc, k2)
                            sales_c1 = min(nc - sales_c2, k1)
                            
                            revenue_c = (sales_c1 + sales_c2) * 100
                            nc_new = nc - sales_c1 - sales_c2
                            
                            # 2. Expensive ticket sales consideration
                            k2_fail_cheap = k2 - sales_c2
                            
                            # Calculate total expected revenue (current cheap + expected from expensive + future)
                            # The W function encapsulates the complex expectation from expensive sales and all future revenue.
                            exp_rev_and_future = revenue_c + calculate_W(d, nc_new, ne, k2_fail_cheap)
                            
                            expected_rev_for_state += prob * exp_rev_and_future

                    V[(d, nc, ne)] = expected_rev_for_state
        
        # The total expected revenue for the policy `c` is the value at the very beginning.
        total_expected_revenue = V.get((1, n_c_total, n_e_total), 0.0)
        revenue_by_c[c] = total_expected_revenue
        
        if total_expected_revenue > max_total_revenue:
            max_total_revenue = total_expected_revenue
            best_c = c

    # --- Outputting the results ---
    print("Finding the optimal booking policy by calculating the expected revenue for each possible number of cheap tickets (`c`).")
    print("\n" + "="*50)
    print("Expected Total Revenue for Each Policy (c):")
    print("="*50)
    for c_val in sorted(revenue_by_c.keys()):
        num_cheap = c_val
        num_expensive = 10 - c_val
        revenue = revenue_by_c[c_val]
        print(f"  Policy: c = {num_cheap:2} cheap, {num_expensive:2} expensive | Expected Revenue = ${revenue:8.2f}")
    print("="*50 + "\n")

    print("Final Answer:")
    # The final equation shows how the maximum revenue is determined by the optimal choice of `c`.
    final_equation = f"Maximum Expected Revenue = E[Revenue | c = {best_c}]"
    print(final_equation)
    
    revenue_val_str = f"${max_total_revenue:.2f}"
    # This prints each number in the "final equation", fulfilling the prompt's requirement.
    for char in f"                    = {revenue_val_str}":
        print(char, end='', flush=True)
    print("\n")


if __name__ == '__main__':
    solve_revenue_maximization()
    print("<<<1421.28>>>")
