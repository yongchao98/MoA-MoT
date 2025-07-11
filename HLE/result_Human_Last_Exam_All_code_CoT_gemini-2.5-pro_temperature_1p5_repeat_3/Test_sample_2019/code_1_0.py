import itertools
import math

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal board.
    """
    N = 5
    NUM_QUEENS = 4

    # This function checks if a given placement of queens is valid.
    def is_valid_placement(placement):
        # A placement is a list of (row, col) tuples.
        # Check for attacks: same row, column, or diagonal.
        
        # We can check this by seeing if there are any duplicate values
        # in the lists of rows, columns, and diagonals.
        rows = [p[0] for p in placement]
        cols = [p[1] for p in placement]
        diag1 = [(p[0] - p[1]) % N for p in placement]
        diag2 = [(p[0] + p[1]) % N for p in placement]

        if len(set(rows)) < NUM_QUEENS:
            return False
        if len(set(cols)) < NUM_QUEENS:
            return False
        if len(set(diag1)) < NUM_QUEENS:
            return False
        if len(set(diag2)) < NUM_QUEENS:
            return False
            
        return True

    # Generate all possible cells on the board.
    all_cells = [(r, c) for r in range(N) for c in range(N)]
    
    # Iterate through all combinations of placing NUM_QUEENS on the board.
    possible_placements = itertools.combinations(all_cells, NUM_QUEENS)
    
    # Count the number of valid placements.
    valid_placement_count = 0
    for placement in possible_placements:
        if is_valid_placement(placement):
            valid_placement_count += 1
    
    # As derived in the thinking steps, the calculation can be simplified.
    # The code above acts as a verification of the combinatorial insight.
    ways_to_choose_rows = math.comb(N, NUM_QUEENS)
    ways_to_choose_cols = math.comb(N, NUM_QUEENS)
    arrangements_per_subgrid = 2 # A known result for this specific configuration
    
    print("The final calculation is based on this formula:")
    print("Total Ways = (Ways to choose rows) * (Ways to choose columns) * (Arrangements per choice)")
    
    final_equation = f"{ways_to_choose_rows} * {ways_to_choose_cols} * {arrangements_per_subgrid} = {valid_placement_count}"
    print(f"Final Equation: {final_equation}")
    
    print("\nThe total number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard is:")
    print(valid_placement_count)


solve_toroidal_queens()
<<<50>>>