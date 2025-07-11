def solve_peaceful_queens():
    """
    Calculates the maximum number m for peaceful white and black queens on a 16x16 board.
    """
    N = 16

    print(f"The board size is {N}x{N}.")
    print("The goal is to find the maximum number m such that m white queens and m black queens can be placed on the board without any white queen attacking a black queen.\n")

    print("Step 1: Prevent diagonal attacks by separating queens by square color.")
    print("We will place all white queens on 'white' squares (where row+col is even).")
    print("We will place all black queens on 'black' squares (where row+col is odd).\n")

    print("Step 2: Prevent row/column attacks by partitioning the board.")
    print(f"We divide the board into quadrants. The size of each quadrant's side is {N} / 2 = {N//2}.")
    print("White queens are restricted to the top-left quadrant (rows 0-7, cols 0-7).")
    print("Black queens are restricted to the bottom-right quadrant (rows 8-15, cols 8-15).\n")

    print("Step 3: Calculate the maximum m.")
    quadrant_side = N // 2
    quadrant_area = quadrant_side * quadrant_side
    
    # In a balanced checkerboard area like an 8x8 quadrant, half the squares are white and half are black.
    available_squares_per_quadrant = quadrant_area // 2

    m = available_squares_per_quadrant

    print("The number of queens 'm' is limited by the number of available squares in each designated zone.")
    print(f"The equation for the number of available squares in a quadrant is: (({N} / 2) * ({N} / 2)) / 2")
    print(f"Calculation: (({quadrant_side}) * ({quadrant_side})) / 2 = {quadrant_area} / 2 = {m}")
    print(f"\nThis means we can place {m} white queens in the top-left quadrant and {m} black queens in the bottom-right quadrant.")
    print(f"Therefore, the maximum value for m is {m}.")

solve_peaceful_queens()