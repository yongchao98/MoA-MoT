def solve_alien_invasion():
    """
    Calculates the maximal size of the alien colony based on a logical deduction.
    """
    board_size = 8 * 8
    print(f"Total squares on the board: {board_size}")

    # To capture the whole board, every stable vacant set must be "punctured"
    # by an initial alien. The smallest stable vacant sets are single rows or columns.
    # This means we must place at least one alien in each row and each column.

    rows_to_cover = 8
    initial_squares = 8
    
    # The two fixed squares d5 and e5 are on the same row.
    rows_covered_by_fixed_squares = 1
    
    # We have 6 other squares to place.
    additional_squares = 6
    
    # How many more rows do we need to cover?
    rows_left_to_cover = rows_to_cover - rows_covered_by_fixed_squares
    print(f"Number of rows that need an initial alien: {rows_left_to_cover}")
    print(f"Number of additional aliens we can place: {additional_squares}")

    # Since we have fewer additional squares than rows that need covering,
    # we cannot cover all rows.
    min_uncovered_rows = rows_left_to_cover - additional_squares
    
    # The smallest stable region we are forced to leave vacant is a single row.
    min_size_of_vacant_set = min_uncovered_rows * 8
    print(f"Minimum size of the final vacant colony: {min_size_of_vacant_set}")
    
    # The maximal colony size is the total board size minus this minimum vacant set.
    max_colony_size = board_size - min_size_of_vacant_set
    
    print("\nThe final calculation is:")
    print(f"K = {board_size} - {min_size_of_vacant_set}")
    print(f"K = {max_colony_size}")
    
solve_alien_invasion()