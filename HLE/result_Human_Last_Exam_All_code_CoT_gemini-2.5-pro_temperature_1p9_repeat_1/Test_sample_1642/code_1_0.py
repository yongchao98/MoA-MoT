def solve_max_queens_problem():
    """
    Solves for the maximum number m of white and black queens on a 16x16 board.
    This is equivalent to solving the 16-Queens problem and dividing the result by 2.
    """
    N = 16

    def is_safe(board, row, col):
        """
        Check if it's safe to place a queen at board[row][col].
        This function is called when 'row' queens have already been placed in
        rows 0 to row-1. We only need to check the left side for attacking queens.
        """
        # board[i] stores the column for the queen in row i.
        for i in range(row):
            # Check for queen in the same column
            if board[i] == col:
                return False
            # Check for queen on the same diagonal
            if abs(board[i] - col) == abs(i - row):
                return False
        return True

    def solve_n_queens_util(board, row):
        """
        A recursive utility function to solve N Queens problem using backtracking.
        """
        # Base case: If all queens are placed, we have a solution.
        if row >= N:
            return True

        # Consider this row and try placing this queen in all columns one by one
        for col in range(N):
            if is_safe(board, row, col):
                # Place this queen in board[row]
                board[row] = col
                
                # Recur to place the rest of the queens
                if solve_n_queens_util(board, row + 1):
                    return True
                
                # If placing queen in board[row] doesn't lead to a solution,
                # then backtrack and remove the queen from board[row]
                # (Overwritten in the next iteration, so no explicit removal needed)

        # If a queen cannot be placed in any column in this row, return false
        return False

    # board[i] will store the column of the queen in row i. Initialize with -1.
    board = [-1] * N
    
    # Note: Finding a solution for N=16 might take a moment.
    # We only need to prove that one solution exists.
    if solve_n_queens_util(board, 0):
        print("This problem is equivalent to finding the maximum number of non-attacking queens on a 16x16 board.")
        print("On an N x N board, the maximum number of non-attacking queens is N.")
        max_queens = N
        print(f"For a 16x16 board, this maximum is {max_queens}.")
        
        print("\nWe are placing m white queens and m black queens.")
        print("The total number of queens is the sum, m + m = 2m.")

        print(f"\nThis gives us the equation: 2 * m = {max_queens}")

        m = max_queens // 2
        print("Solving for m, we divide the total number of queens by 2:")
        print(f"m = {max_queens} / 2 = {m}")
        
        print(f"\nThus, the maximum number m is {m}.")

    else:
        # This part should not be reached for N>=4
        print("It was not possible to place 16 queens on the board.")

if __name__ == "__main__":
    solve_max_queens_problem()