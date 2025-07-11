import sys
from functools import lru_cache

def solve_trading_problem(initial_investment, num_trades):
    """
    Calculates the maximum probability of doubling an initial investment
    in a specific number of trades using dynamic programming.

    Args:
        initial_investment (int): The starting amount of money (M).
        num_trades (int): The exact number of trades to be made (n).
    """
    # We need to increase the recursion limit for problems with a large number of trades.
    # The default limit might be too low.
    if num_trades + 100 > sys.getrecursionlimit():
        sys.setrecursionlimit(num_trades + 100)
    
    TARGET_INVESTMENT = 2 * initial_investment

    # Use a Least Recently Used (LRU) cache for memoization.
    # This decorator automatically stores the results of the function calls,
    # preventing re-computation for the same (m, l) state.
    @lru_cache(maxsize=None)
    def find_max_prob(money, trades_left):
        """
        The core recursive function to solve the subproblem.

        Args:
            money (int): The current amount of money.
            trades_left (int): The number of trades remaining.

        Returns:
            float: The maximum probability of success from this state.
        """
        # Base Case: No trades left.
        # Success if money is exactly the target, failure otherwise.
        if trades_left == 0:
            return 1.0 if money == TARGET_INVESTMENT else 0.0

        # Base Case: No money to make a trade (or capital wiped out).
        # We can't proceed, so probability of reaching the target is 0.
        if money <= 0:
            return 0.0

        # --- Calculate probability of success if choosing Strategy Alpha ---
        prob_alpha = 0.0
        fee_alpha = 1
        # Check if we can afford the investment fee for Strategy Alpha.
        if money >= fee_alpha:
            # The probability equation for taking this strategy path.
            # 60% chance of success, returning 2 (net profit 1).
            # 40% chance of failure, returning 0 (net loss 1).
            prob_alpha = (0.60 * find_max_prob(money + 1, trades_left - 1) +
                          0.40 * find_max_prob(money - 1, trades_left - 1))

        # --- Calculate probability of success if choosing Strategy Beta ---
        prob_beta = 0.0
        fee_beta = 3
        # Check if we can afford the investment fee for Strategy Beta.
        if money >= fee_beta:
            # The probability equation for taking this strategy path.
            # 20% chance of success, returning 15 (net profit 12).
            # 80% chance of failure, returning 0 (net loss 3).
            prob_beta = (0.20 * find_max_prob(money + 12, trades_left - 1) +
                         0.80 * find_max_prob(money - 3, trades_left - 1))

        # The optimal strategy at state (money, trades_left) is to choose
        # the action (Alpha or Beta) that yields the highest probability of success.
        return max(prob_alpha, prob_beta)

    # Initial call to the recursive function with starting parameters.
    max_probability = find_max_prob(initial_investment, num_trades)
    
    print(f"For an initial investment M = {initial_investment} and n = {num_trades} trades:")
    print(f"The maximum probability of achieving exactly {TARGET_INVESTMENT} GBP is: {max_probability}")


if __name__ == '__main__':
    # Define the inputs for the problem
    M = 25  # Initial Investment
    n = 10   # Number of Trades
    solve_trading_problem(M, n)