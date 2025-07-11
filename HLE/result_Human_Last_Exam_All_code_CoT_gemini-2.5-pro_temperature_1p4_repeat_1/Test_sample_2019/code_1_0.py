import sys

# It is recommended to use PyPy for this script for better performance,
# though standard Python will also work.

def solve():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4

    # memoization cache for the recursive function
    memo = {}

    def is_safe(r, c, placements):
        """
        Checks if a queen at (r, c) attacks any queen in the placements list.
        """
        for pr, pc in placements:
            if r == pr:
                return False
            if c == pc:
                return False
            # Check for toroidal diagonal attacks
            if (r - c) % N == (pr - pc) % N:
                return False
            if (r + c) % N == (pr + pc) % N:
                return False
        return True

    def count_placements(k, start_index, placements):
        """
        Recursively counts valid queen placements using backtracking.
        k: number of queens left to place.
        start_index: board square index to start searching from (avoids duplicates).
        placements: a tuple of (r, c) for already placed queens.
        """
        # Base case: successfully placed all queens
        if k == 0:
            return 1
        
        # Memoization key
        state = (k, start_index, placements)
        if state in memo:
            return memo[state]

        count = 0
        # Iterate through all possible squares starting from start_index
        for i in range(start_index, N * N):
            # Optimization: if remaining squares are fewer than queens to place, stop.
            if (N * N - i) < k:
                break
                
            r = i // N
            c = i % N

            if is_safe(r, c, placements):
                # If safe, place queen and recurse for the next one
                new_placements = placements + ((r, c),)
                count += count_placements(k - 1, i + 1, new_placements)
        
        memo[state] = count
        return count

    # Initial call to the recursive function
    total_ways = count_placements(NUM_QUEENS, 0, tuple())

    # Formulate the output as an equation based on symmetry, as explained in the plan.
    # C = (Total * Num_Queens) / Board_Area
    # This C represents the number of ways to place the other 3 queens
    # once one has been fixed on a single square.
    
    board_area = N * N
    c = (total_ways * NUM_QUEENS) / board_area

    if c == int(c):
      c = int(c)
      
    print(f"Based on the symmetry of the toroidal board, we can establish the following relationship:")
    print(f"({c} * {board_area}) / {NUM_QUEENS} = {total_ways}")
    print(f"\nThis means fixing one queen on any square leaves {c} ways to place the other 3.")
    print(f"The program confirms the total number of ways is {total_ways}.")


solve()
<<<50>>>