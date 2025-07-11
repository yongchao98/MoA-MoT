import numpy as np

def solve_airline_revenue():
    """
    Solves the airline revenue management problem using Dynamic Programming.
    """
    
    # DP table: V[d][s] = max expected revenue with d days left and s seats.
    # d ranges from 0 to 14. s ranges from 0 to 10.
    V = np.zeros((15, 11))
    
    # policy[d][s] stores the optimal protection level k for state (d, s)
    policy = np.zeros((15, 11), dtype=int)

    # Probability distributions for customer arrivals
    P1 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 1 demands
    P2 = {0: 0.25, 1: 0.5, 2: 0.25}  # Class 2 demands

    # Recursive helper function to calculate the outcome of a single day's sales
    def solve_day(seats, n1, n2, k, V_prev_day, memo):
        """
        Calculates the expected day revenue and expected future revenue for a
        given state and customer arrivals, handling sequential, probabilistic sales.
        
        Args:
            seats (int): Number of seats available at the start of the day.
            n1 (int): Number of Class 1 customer arrivals.
            n2 (int): Number of Class 2 customer arrivals.
            k (int): Protection level for the day.
            V_prev_day (np.array): DP values for the next day (d-1).
            memo (dict): Memoization table for the recursive calls.

        Returns:
            tuple: (expected_day_revenue, expected_future_revenue)
        """
        state = (seats, n1, n2)
        if state in memo:
            return memo[state]

        # Base case: no more customers or seats. Return current state's future value.
        if (n1 == 0 and n2 == 0) or seats == 0:
            return (0, V_prev_day[seats])

        # --- Process customers sequentially: Class 2 has priority ---
        if n2 > 0:
            if seats > k:  # Cheap ticket is available
                day_rev_cont, v_next_cont = solve_day(seats - 1, n1, n2 - 1, k, V_prev_day, memo)
                result = (100 + day_rev_cont, v_next_cont)
            else:  # Offer expensive ticket (if seats > 0, which is true here)
                # Case A: Customer BUYS (prob 0.5)
                rev_A, v_next_A = solve_day(seats - 1, n1, n2 - 1, k, V_prev_day, memo)
                
                # Case B: Customer WALKS AWAY (prob 0.5)
                rev_B, v_next_B = solve_day(seats, n1, n2 - 1, k, V_prev_day, memo)
                
                exp_day_rev = 0.5 * (200 + rev_A) + 0.5 * rev_B
                exp_v_next = 0.5 * v_next_A + 0.5 * v_next_B
                result = (exp_day_rev, exp_v_next)
        else:  # Process Class 1 customers (n2 is 0, n1 > 0)
            if seats > k:  # Cheap ticket is available
                day_rev_cont, v_next_cont = solve_day(seats - 1, n1 - 1, 0, k, V_prev_day, memo)
                result = (100 + day_rev_cont, v_next_cont)
            else:  # C1 customer walks away
                day_rev_cont, v_next_cont = solve_day(seats, n1 - 1, 0, k, V_prev_day, memo)
                result = (day_rev_cont, v_next_cont)
        
        memo[state] = result
        return result

    # --- Main DP loop, iterating backwards from the last day ---
    for d in range(1, 15):  # d = days remaining (1 to 14)
        for s in range(1, 11):  # s = seats available (1 to 10)
            max_exp_revenue = -1
            best_k = -1
            
            # Find the best protection level k for this state (d, s)
            for k in range(s + 1):
                memo_solver = {}
                current_k_total_exp_rev = 0
                
                # Determine customer arrival distributions for this day
                # Last 7 days (d=1 to 7) are week 2
                is_week_2 = (d <= 7)
                d2_dist = P2 if is_week_2 else {0: 1.0}
                
                # Calculate expected value by summing over all arrival scenarios
                for d1, p1 in P1.items():
                    for d2, p2 in d2_dist.items():
                        prob = p1 * p2
                        day_rev, future_rev = solve_day(s, d1, d2, k, V[d - 1], memo_solver)
                        current_k_total_exp_rev += prob * (day_rev + future_rev)

                if current_k_total_exp_rev > max_exp_revenue:
                    max_exp_revenue = current_k_total_exp_rev
                    best_k = k
            
            V[d][s] = max_exp_revenue
            policy[d][s] = best_k

    # --- Final Result Calculation and Output ---
    final_revenue = V[14][10]
    k_opt = policy[14][10]

    print("Finding the maximum expected total revenue for 10 seats over 14 days.")
    print(f"The optimal policy on the first day (Day 14) with 10 seats is to protect k={k_opt} seats for expensive fares.\n")
    print("The maximum expected revenue is calculated as the sum of outcomes for all possible customer demands on Day 14:")
    print("V(14, 10) = Sum[ P(d1) * (E[Day 14 Revenue | d1] + E[Future Revenue V(13) | d1]) ]\n")

    total_check = 0
    # Day 14 is in the first week, so only Class 1 customers arrive (d2=0)
    for d1, p1 in P1.items():
        memo_final = {}
        day_rev, future_rev = solve_day(10, d1, 0, k_opt, V[13], memo_final)
        term_value = day_rev + future_rev
        total_check += p1 * term_value
        print(f"Term for d1={d1} (Prob={p1:.2f}):")
        print(f"  {p1:.2f} * ({day_rev:9.2f} + {future_rev:9.2f}) = {p1 * term_value:9.2f}")

    print("-" * 50)
    print(f"Total Maximum Expected Revenue = {total_check:.2f}\n")
    print("Final Answer:")
    print(f"{final_revenue:.2f}")
    
    return final_revenue

# Run the calculation and get the final answer
final_answer = solve_airline_revenue()
# Present the final numerical answer in the required format
print(f"<<<{final_answer:.2f}>>>")
