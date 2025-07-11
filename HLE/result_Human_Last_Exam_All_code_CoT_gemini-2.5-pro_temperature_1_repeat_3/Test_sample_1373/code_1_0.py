import sys

# Since the number of trades 'n' can be large, the recursion depth might
# exceed Python's default limit. We increase it to be safe.
sys.setrecursionlimit(3000)

# Use a dictionary for memoization to store results of subproblems.
# The key will be a tuple (m, l) representing the state.
memo = {}

# These will be set by the main function call.
INITIAL_MONEY = 0
TOTAL_TRADES = 0

def find_max_probability(current_money, trades_done):
    """
    Calculates the maximum probability of reaching the target investment using dynamic programming.

    Args:
        current_money (int): The amount of money currently held.
        trades_done (int): The number of trades already executed.

    Returns:
        float: The maximum probability of success from this state.
    """
    global INITIAL_MONEY, TOTAL_TRADES, memo

    # Base Case: If all trades are completed.
    if trades_done == TOTAL_TRADES:
        # Success if money is exactly double the initial investment, otherwise failure.
        return 1.0 if current_money == 2 * INITIAL_MONEY else 0.0

    # Memoization: If we have already solved this subproblem, return the stored result.
    state = (current_money, trades_done)
    if state in memo:
        return memo[state]

    # --- Recursive Step: Calculate probabilities for each possible strategy ---

    # Initialize probabilities for the two strategies. -1.0 indicates not possible.
    prob_if_alpha = -1.0
    prob_if_beta = -1.0

    # Strategy Alpha: Requires £1 fee.
    if current_money >= 1:
        # 60% chance of success (return £2 -> net profit £1)
        # 40% chance of failure (return £0 -> net loss £1)
        prob_if_alpha = (0.6 * find_max_probability(current_money + 1, trades_done + 1) +
                         0.4 * find_max_probability(current_money - 1, trades_done + 1))

    # Strategy Beta: Requires £3 fee.
    if current_money >= 3:
        # 20% chance of success (return £15 -> net profit £12)
        # 80% chance of failure (return £0 -> net loss £3)
        prob_if_beta = (0.2 * find_max_probability(current_money + 12, trades_done + 1) +
                        0.8 * find_max_probability(current_money - 3, trades_done + 1))

    # The optimal strategy is the one that maximizes the probability of success.
    # If a strategy is not affordable, its probability is -1.0, so max() will pick the valid one.
    # If neither is affordable, max will be -1.0, which we correct to 0.0, as it's a failed path.
    result = max(prob_if_alpha, prob_if_beta)
    if result < 0:
        result = 0.0
        
    # Store the result in the memoization table before returning.
    memo[state] = result
    return result

def solve_trading_problem(M, n):
    """
    A wrapper function to set global parameters, run the calculation, and print the result.
    """
    global INITIAL_MONEY, TOTAL_TRADES, memo
    INITIAL_MONEY = M
    TOTAL_TRADES = n
    # Clear the memoization table for a fresh run.
    memo.clear()

    # The problem starts with M pounds and 0 trades done.
    final_probability = find_max_probability(M, 0)
    
    target_money = 2 * M

    print(f"Initial investment M: {M} GBP")
    print(f"Number of trades n: {n}")
    print(f"Target investment: {target_money} GBP")
    print(f"The maximum probability of reaching exactly {target_money} GBP after {n} trades is: {final_probability}")


# --- Main execution with example values ---
if __name__ == '__main__':
    # You can change these values to test different scenarios.
    M_val = 25
    n_val = 10
    solve_trading_problem(M_val, n_val)