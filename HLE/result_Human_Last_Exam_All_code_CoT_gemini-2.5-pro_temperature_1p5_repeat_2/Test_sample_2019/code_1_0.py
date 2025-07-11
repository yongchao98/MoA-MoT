import itertools
from math import comb

def count_valid_permutations(rows, cols, board_size=5):
    """
    Counts the number of valid queen placements for a fixed set of rows and columns.
    A placement is valid if no two queens attack each other on the diagonals.
    
    Args:
        rows (tuple): A tuple of row indices.
        cols (tuple): A tuple of column indices.
        board_size (int): The size of the board (N for an NxN board).

    Returns:
        int: The number of valid permutations.
    """
    num_queens = len(rows)
    valid_perms_count = 0
    
    # Iterate through all possible assignments of columns to the fixed rows
    for p_cols in itertools.permutations(cols):
        # A placement is a list of (row, column) coordinates
        placement = list(zip(rows, p_cols))
        
        # Check for main diagonal attacks
        main_diags = set((r - c) % board_size for r, c in placement)
        if len(main_diags) < num_queens:
            continue
            
        # Check for anti-diagonal attacks
        anti_diags = set((r + c) % board_size for r, c in placement)
        if len(anti_diags) < num_queens:
            continue
            
        # If no diagonal attacks, this is a valid permutation
        valid_perms_count += 1
        
    return valid_perms_count

# Define board parameters
BOARD_SIZE = 5
NUM_QUEENS = 4

# 1. Calculate the number of ways to choose the rows and columns
num_row_choices = comb(BOARD_SIZE, NUM_QUEENS)
num_col_choices = comb(BOARD_SIZE, NUM_QUEENS)

# 2. Calculate the number of valid permutations for a canonical choice of rows/cols
# We can use rows {0, 1, 2, 3} and columns {0, 1, 2, 3} as a representative case.
canonical_rows = tuple(range(NUM_QUEENS))
canonical_cols = tuple(range(NUM_QUEENS))
num_perms = count_valid_permutations(canonical_rows, canonical_cols, BOARD_SIZE)

# 3. Calculate the total number of ways
total_ways = num_row_choices * num_col_choices * num_perms

# Print the breakdown of the calculation as requested
print("The calculation for the total number of ways is as follows:")
print(f"Ways to choose {NUM_QUEENS} rows from {BOARD_SIZE}: {num_row_choices}")
print(f"Ways to choose {NUM_QUEENS} columns from {BOARD_SIZE}: {num_col_choices}")
print(f"Valid permutations for any choice of {NUM_QUEENS} rows and columns: {num_perms}")
print("\nFinal calculation:")
print(f"{num_row_choices} * {num_col_choices} * {num_perms} = {total_ways}")