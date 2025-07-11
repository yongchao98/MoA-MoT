def solve_max_queens():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 chessboard without attacking each other.
    The script explains the logical steps to arrive at the solution.
    """
    # Define the size of the chessboard
    N = 16

    print(f"The problem is set on a {N}x{N} chessboard.")
    print("We need to find the maximum number m for two sets of queens (m white, m black).")
    print("The condition 'without attacking each other' means no queen can attack any other queen on the board.")
    
    print("\nStep 1: Calculate the total number of queens.")
    print("Total queens = (m white queens) + (m black queens) = 2 * m.")
    
    print("\nStep 2: Apply the N-Queens problem principle.")
    print(f"A key rule in chess problems is that on an {N}x{N} board, you can place at most {N} non-attacking queens.")
    print(f"Therefore, the total number of queens (2 * m) must be less than or equal to {N}.")
    
    print("\nStep 3: Formulate and solve the inequality.")
    print("This gives us the inequality including the final numbers:")
    print(f"2 * m <= {N}")
    
    print("\nSolving for m:")
    max_m = N // 2
    print(f"m <= {N} / 2")
    print(f"m <= {max_m}")
    
    print(f"\nConclusion: The maximum number m is {max_m}.")

solve_max_queens()