import sys

# It's a recursive problem, so we might hit the recursion limit for large n.
# Increase the recursion limit for safety.
sys.setrecursionlimit(2000)

# Memoization table (dictionary) to store results of subproblems
memo = {}
TARGET_MONEY = 0

def solve(current_money, trades_left):
    """
    Calculates the maximum probability of success using dynamic programming.

    Args:
        current_money (int): The amount of money at the current state.
        trades_left (int): The number of trades remaining.

    Returns:
        float: The maximum probability of reaching TARGET_MONEY.
    """
    global TARGET_MONEY, memo

    # Base Case: No trades left
    if trades_left == 0:
        return 1.0 if current_money == TARGET_MONEY else 0.0

    # Memoization Check: If we have already solved this subproblem, return the stored result
    if (current_money, trades_left) in memo:
        return memo[(current_money, trades_left)]
        
    # If money is not enough to make any trade, the probability of reaching the target is 0
    # (unless we are in the base case, which is handled above).
    if current_money < 1:
        return 0.0

    # --- Calculate probability for Strategy Alpha ---
    prob_alpha = 0.0
    # Strategy Alpha is only possible if we have at least £1
    if current_money >= 1:
        # Prob_success = P(win) * solve(state_after_win) + P(loss) * solve(state_after_loss)
        prob_alpha = (0.6 * solve(current_money + 1, trades_left - 1) +
                      0.4 * solve(current_money - 1, trades_left - 1))

    # --- Calculate probability for Strategy Beta ---
    prob_beta = 0.0
    # Strategy Beta is only possible if we have at least £3
    if current_money >= 3:
        prob_beta = (0.2 * solve(current_money + 12, trades_left - 1) +
                     0.8 * solve(current_money - 3, trades_left - 1))

    # The optimal choice is the one with the maximum probability of success
    result = max(prob_alpha, prob_beta)

    # Store the result in the memoization table before returning
    memo[(current_money, trades_left)] = result
    
    return result

def main():
    """
    Main function to set parameters and run the simulation.
    """
    global TARGET_MONEY, memo
    
    # --- Parameters ---
    # Initial Investment (in GBP)
    M = 25
    # Number of trades to be executed
    n = 10

    # The goal is to exactly double the initial investment
    TARGET_MONEY = 2 * M
    
    # Clear the memoization table for a fresh run
    memo.clear()

    # Calculate the overall maximum probability starting with M money and n trades
    # Note: The 'solve' function populates the memo table as a side effect
    final_probability = solve(M, n)

    # --- Output the results and the "final equation" ---
    print(f"Initial State: Money = £{M}, Trades = {n}, Target = £{TARGET_MONEY}\n")

    # To show the final equation, we determine the optimal first move and its components
    if M < 1:
        print("No trades are possible with the initial investment.")
        print(f"P({M}, {n}) = 0.0")
        return

    # Probability if the first move is Alpha
    # The required subproblem results are already in the memo table
    prob_alpha_first_step = (0.6 * solve(M + 1, n - 1) +
                             0.4 * solve(M - 1, n - 1))

    # Determine which strategy is optimal for the first trade
    if M < 3: # Beta is not affordable, Alpha is the only choice
        print("Optimal first move is Strategy Alpha (only affordable option).")
        p_win = memo.get((M + 1, n - 1), 0)
        p_loss = memo.get((M - 1, n - 1), 0)
        print(f"P({M}, {n}) = 0.60 * P({M + 1}, {n - 1}) + 0.40 * P({M - 1}, {n - 1})")
        print(f"P({M}, {n}) = 0.60 * {p_win:.5f} + 0.40 * {p_loss:.5f}")
    else: # Both strategies are affordable, we must compare them
        prob_beta_first_step = (0.2 * solve(M + 12, n - 1) +
                                0.8 * solve(M - 3, n - 1))
        
        if prob_alpha_first_step >= prob_beta_first_step:
            print("Optimal first move is Strategy Alpha.")
            p_win = memo.get((M + 1, n - 1), 0)
            p_loss = memo.get((M - 1, n - 1), 0)
            print(f"P({M}, {n}) = 0.60 * P({M + 1}, {n - 1}) + 0.40 * P({M - 1}, {n - 1})")
            print(f"P({M}, {n}) = 0.60 * {p_win:.5f} + 0.40 * {p_loss:.5f}")
        else:
            print("Optimal first move is Strategy Beta.")
            p_win = memo.get((M + 12, n - 1), 0)
            p_loss = memo.get((M - 3, n - 1), 0)
            print(f"P({M}, {n}) = 0.20 * P({M + 12}, {n - 1}) + 0.80 * P({M - 3}, {n - 1})")
            print(f"P({M}, {n}) = 0.20 * {p_win:.5f} + 0.80 * {p_loss:.5f}")

    print("-" * 30)
    print(f"Maximum probability of success: {final_probability:.5f}")
    print("-" * 30)


if __name__ == '__main__':
    main()