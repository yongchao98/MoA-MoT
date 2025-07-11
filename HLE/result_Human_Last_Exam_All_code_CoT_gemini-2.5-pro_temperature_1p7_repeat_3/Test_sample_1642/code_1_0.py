def solve_queen_problem():
    """
    This script calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    
    # Board dimension
    N = 16
    
    print("The problem is to find the maximum number 'm' for 'm' white and 'm' black queens on a 16x16 board.")
    print("The condition is that no queen attacks any other queen, regardless of color.")
    print("\nThis simplifies the problem to finding the maximum total number of non-attacking queens on a 16x16 board.")
    
    # On an N x N board, the maximum number of non-attacking queens is N.
    max_queens = N
    
    print(f"\nOn a {N}x{N} board, you can place at most {max_queens} non-attacking queens, as each requires a unique row.")
    
    # The total number of queens is m + m = 2m
    print("\nThe total number of queens is the sum of white (m) and black (m) queens, which is 2 * m.")
    
    print(f"\nWe can set up the equation: 2 * m = (maximum number of queens)")
    
    # Define the numbers in the final equation
    coefficient_m = 2
    total_queens = max_queens
    
    print("\nThe final equation is:")
    print(f"{coefficient_m} * m = {total_queens}")
    
    # Solve for m
    m = total_queens // coefficient_m
    
    print("\nSolving for m:")
    print(f"m = {total_queens} / {coefficient_m}")
    print(f"m = {m}")
    
    print(f"\nTherefore, the maximum number m is {m}.")

solve_queen_problem()
<<<8>>>