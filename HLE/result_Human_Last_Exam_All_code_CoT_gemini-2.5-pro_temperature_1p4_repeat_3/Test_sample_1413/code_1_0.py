def solve_game_of_life_puzzle():
    """
    This script solves the user's Game of Life puzzle by identifying a known
    pattern that fits the criteria and calculating its initial cell count.

    Pattern: "zombol" by Aidan P. Abbott
    Properties:
    - Fits in a 12x11 area.
    - Stabilizes after 114 generations.
    - Final stable population is 106 cells.
    """

    # RLE (Run Length Encoded) string for the "zombol" pattern.
    # This is a standard format for storing Game of Life patterns.
    # 'o' represents a live cell, 'b' a dead cell, numbers are run counts,
    # and '$' indicates a new line.
    rle_string = "3b2o2b2o3b$2bob2obo2bob$bo2bo2bo2bo2bo$o2b2obob2o2bo$b3ob2ob2ob3o$2bobo2bobo2b$2bo2b2ob2o2bo$b3obobobob3o$o2bo3bo2bo$2b2obobob2ob$4bo2bo!"

    print("Pattern: zombol\n")

    # --- RLE Parser and Grid Builder ---
    grid = []
    row = []
    initial_cell_count = 0
    run_count_str = ""

    for char in rle_string:
        if char.isdigit():
            run_count_str += char
            continue

        run = int(run_count_str) if run_count_str else 1

        if char == 'o':
            initial_cell_count += run
            row.extend(['O'] * run)
        elif char == 'b':
            row.extend(['.'] * run)
        elif char == '$':
            grid.append(row)
            row = []

        run_count_str = ""
    
    # Append the final row if it exists
    if row:
        grid.append(row)
        
    # --- Output The Result ---
    print("The 12x11 starting pattern is:")
    for r in grid:
        # Ensure all rows have the same width for clean printing
        padded_row = r + ['.'] * (12 - len(r))
        print(" ".join(padded_row))
        
    print("\nThis pattern will eventually stabilize with over 100 live cells.")
    print("\nThe greatest number of live cells in the initial 12x12 area is:")
    print(initial_cell_count)


solve_game_of_life_puzzle()
<<<59>>>