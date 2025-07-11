def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    """
    
    # Size of the chessboard
    N = 16
    
    print(f"The user wants to find the maximum number m for a {N}x{N} board.")
    print("The condition 'without attacking each other' means no queen can attack any other queen, regardless of color.")
    print(f"This means we are placing a total of m + m = 2*m queens on a {N}x{N} board.")
    print("\nStep 1: State the maximum number of non-attacking queens on an NxN board.")
    print(f"For an N x N board (where N >= 4), the maximum number of non-attacking queens is N.")
    print(f"For our {N}x{N} board, the maximum total number of queens is {N}.")
    
    # Total number of queens
    total_queens = N
    
    print("\nStep 2: Set up the equation.")
    print(f"The total number of queens is 2 * m. This total cannot exceed the maximum of {N}.")
    print(f"So, the equation is: 2 * m = {total_queens}")

    print("\nStep 3: Solve for m.")
    # In Python, we use floating point division by default.
    # To ensure integer arithmetic is clear, we can use //.
    m = total_queens // 2
    
    print(f"m = {total_queens} / 2")
    print(f"m = {m}")
    
    print(f"\nThus, the maximum number m is {m}.")

solve_queen_problem()