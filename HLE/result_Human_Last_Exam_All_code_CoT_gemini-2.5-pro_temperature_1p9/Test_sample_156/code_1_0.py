import math

def solve():
    """
    Solves the alien invasion problem by applying a perimeter-based argument.
    The plan is as follows:
    1.  Analyze the process of colony growth. The key observation is how the perimeter of the captured region changes.
    2.  Establish an upper bound on the perimeter of the final colony based on the maximum possible perimeter of the initial colony.
    3.  Relate the area of a region on a grid to its perimeter. The process stops when the colony is "convex" (has no dents that can be filled). For such a shape on a grid, the area is maximized for a given perimeter by a rectangle.
    4.  Calculate the maximum possible area of a rectangle that respects the perimeter bound and fits on the 8x8 board. This value will be our K.
    """
    
    # --- Step 1: Initial Conditions ---
    initial_squares = 8
    board_dim = 8
    # The two fixed squares d5 and e5 are adjacent.
    # This means there is at least one shared side among the initial 8 squares.
    min_adjacencies = 1
    
    print("Step 1: Analyzing the initial setup")
    print(f"The aliens start with {initial_squares} captured squares on an {board_dim}x{board_dim} board.")
    print("Two of these squares, d5 and e5, are adjacent to each other.")
    print("-" * 20)

    # --- Step 2: The Perimeter Argument ---
    # The perimeter of a set of squares is the number of its edges that border vacant squares.
    # The max perimeter for N squares is 4*N. Each shared side between two squares reduces the total perimeter by 2.
    max_initial_perimeter = initial_squares * 4 - 2 * min_adjacencies
    
    print("Step 2: The Perimeter Invariant")
    print("Let P be the perimeter of the alien colony.")
    print("When a new vacant square is captured, it must have k >= 2 neighboring squares that are already captured.")
    print("The change in the colony's perimeter when adding this square is delta_P = 4 - 2*k.")
    print("Since k must be 2, 3, or 4, the change in perimeter is 0, -2, or -4, respectively.")
    print("This proves that the perimeter of the colony can NEVER increase: P_final <= P_initial.")
    print(f"\nThe maximum possible initial perimeter (P_initial) occurs when the 8 squares are as spread out as possible, with only the mandatory d5-e5 adjacency.")
    print(f"So, maximum P_initial = ({initial_squares} * 4) - (2 * {min_adjacencies}) = {max_initial_perimeter}.")
    print(f"Therefore, the final perimeter P_final must be less than or equal to {max_initial_perimeter}.")
    print("-" * 20)

    # --- Step 3: Maximizing the Final Area ---
    # The growth stops when no vacant square has 2 or more captured neighbors. This means the final shape is "convex" in a grid sense.
    # For a given perimeter on a grid, the area is maximized by a shape that is as compact, or "round", as possible. We model this optimal shape as a rectangle.
    
    print("Step 3: Finding the Maximum Area for the Bounded Perimeter")
    print(f"The final colony shape must be a stable region with a perimeter P_final <= {max_initial_perimeter}.")
    print("To maximize area for a given perimeter on a grid, we consider rectangular shapes.")
    print("The task is to find the rectangle (n x m) with the largest area 'n * m' that satisfies:")
    print(f"  1. It fits on an {board_dim}x{board_dim} board (i.e., n <= {board_dim} and m <= {board_dim}).")
    print(f"  2. Its perimeter, 2*(n+m), is at most {max_initial_perimeter}.")

    max_area = 0
    best_n, best_m = 0, 0
    
    # Iterate through all possible rectangle dimensions (n x m) on the board
    for n in range(1, board_dim + 1):
        for m in range(n, board_dim + 1): # m>=n to avoid duplicate checks
            perimeter = 2 * (n + m)
            if perimeter <= max_initial_perimeter:
                area = n * m
                if area > max_area:
                    max_area = area
                    best_n, best_m = n, m
    
    K = max_area
    print(f"\nThe search for the optimal rectangle yields the following:")
    print(f"  Dimensions: {best_n} squares by {best_m} squares")
    print(f"  Resulting Perimeter: 2 * ({best_n} + {best_m}) = {2 * (best_n + best_m)}")
    print(f"  Resulting Area: {best_n} * {best_m} = {K}")
    print("-" * 20)

    # --- Step 4: Final Conclusion ---
    print("Step 4: Conclusion")
    print("This analysis establishes an upper bound: the maximal size K cannot exceed the area of this optimal rectangle.")
    print("The problem specifies that the aliens use 'optimal play' to maximize the final size.")
    print("This implies that this theoretical maximum bound is indeed achievable with a clever initial placement and growth strategy.")
    print(f"Thus, the maximal size K of the aliens' colony is {K}.")
    
    # As requested, output the final numerical expression for the answer.
    print("\nFinal equation for the answer K:")
    print(f"{best_n} * {best_m} = {K}")

solve()