import sys

def solve_trading_problem(initial_investment_M, num_trades_n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades.

    Args:
        initial_investment_M (int): The starting capital.
        num_trades_n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of success.
    """
    # Increase recursion limit for larger n, as the problem can have long dependency chains.
    # This is not strictly needed for the iterative DP but good practice for recursive versions.
    # sys.setrecursionlimit(max(sys.getrecursionlimit(), 2 * num_trades_n))

    # --- Constants from the problem description ---
    # Strategy Alpha
    FEE_ALPHA = 1
    NET_PROFIT_ALPHA = 1 # Net profit: +1 on success, -1 on failure
    PROB_SUCCESS_ALPHA = 0.6
    
    # Strategy Beta
    FEE_BETA = 3
    NET_PROFIT_BETA = 12 # Net profit: +12 on success, -3 on failure
    PROB_SUCCESS_BETA = 0.2
    
    # --- DP setup ---
    target_capital = 2 * initial_investment_M
    
    # A safe upper bound for capital to consider. If capital is m with k trades
    # left, to reach 2M we need 2M <= m + k*12, so m >= 2M - 12k.
    # Also, we need 2M >= m - k*3, so m <= 2M + 3k.
    # The max value of m occurs for max k=n, so m <= 2M + 3n.
    max_capital = target_capital + 3 * num_trades_n
    
    # DP table (space-optimized). dp_next corresponds to probabilities for l+1 trades.
    dp_next = [0.0] * (max_capital + 1)
    
    # Base case: l = n (after n trades are complete)
    # Success only if capital is exactly target_capital.
    if target_capital <= max_capital:
        dp_next[target_capital] = 1.0
        
    # Iterate backwards from l = n-1 down to 0
    for l in range(num_trades_n - 1, -1, -1):
        dp_curr = [0.0] * (max_capital + 1)
        for m in range(max_capital + 1):
            
            prob_alpha = 0.0
            if m >= FEE_ALPHA:
                # Calculate probability of success if choosing Alpha
                prob_succ_outcome = 0.0
                if m + NET_PROFIT_ALPHA <= max_capital:
                    prob_succ_outcome = dp_next[m + NET_PROFIT_ALPHA]
                
                prob_fail_outcome = 0.0
                if m - FEE_ALPHA >= 0:
                    prob_fail_outcome = dp_next[m - FEE_ALPHA]
                
                prob_alpha = PROB_SUCCESS_ALPHA * prob_succ_outcome + (1 - PROB_SUCCESS_ALPHA) * prob_fail_outcome
            
            prob_beta = 0.0
            if m >= FEE_BETA:
                # Calculate probability of success if choosing Beta
                prob_succ_outcome = 0.0
                if m + NET_PROFIT_BETA <= max_capital:
                    prob_succ_outcome = dp_next[m + NET_PROFIT_BETA]

                prob_fail_outcome = 0.0
                if m - FEE_BETA >= 0:
                    prob_fail_outcome = dp_next[m - FEE_BETA]
                
                prob_beta = PROB_SUCCESS_BETA * prob_succ_outcome + (1 - PROB_SUCCESS_BETA) * prob_fail_outcome

            # Choose the optimal strategy
            dp_curr[m] = max(prob_alpha, prob_beta)
            
        # Current row becomes next row for the previous iteration
        dp_next = dp_curr

    # The final answer is the probability at the start (l=0, m=M)
    result = dp_next[initial_investment_M]
    return result

if __name__ == '__main__':
    # Example values for M and n
    M = 25
    n = 10
    
    probability = solve_trading_problem(M, n)
    
    print(f"Starting with M = {M} and making n = {n} trades.")
    print(f"The maximum probability of reaching exactly {2*M} is: {probability}")

    # Another example from the thought process
    M_test = 1
    n_test = 1
    prob_test = solve_trading_problem(M_test, n_test)
    print(f"\nStarting with M = {M_test} and making n = {n_test} trades.")
    # Expected: use Alpha for 0.6 probability of getting m+1=2.
    print(f"The maximum probability of reaching exactly {2*M_test} is: {prob_test}")

