def solve_queen_problem():
    """
    Calculates the maximum number m such that m white queens and m black queens
    can coexist on a 16x16 board without attacking each other.
    """
    # Define the size of the chessboard
    n = 16

    print(f"The chessboard size is {n}x{n}.")
    print("We need to place 'm' white queens and 'm' black queens.")
    print("This makes for a total of 2 * m queens on the board.")
    print("\nThe core rule is that no queen can attack any other queen.")
    print(f"This is the classic N-Queens problem, where N={n}.")
    print(f"\nFor an {n}x{n} board, the maximum number of non-attacking queens is {n}.")
    print("So, the total number of queens (2 * m) must be equal to this maximum number.")
    
    # Formulate and display the equation
    print("\nThis gives us the final equation:")
    # The variable for the total number of queens is 2*m
    # The maximum number of non-attacking queens on an n*n board is n
    print(f"2 * m = {n}")

    # Solve for m
    m = n // 2
    
    print("\nSolving for m:")
    print(f"m = {n} / 2")
    print(f"m = {m}")

    print(f"\nTherefore, the maximum number m is {m}.")

solve_queen_problem()