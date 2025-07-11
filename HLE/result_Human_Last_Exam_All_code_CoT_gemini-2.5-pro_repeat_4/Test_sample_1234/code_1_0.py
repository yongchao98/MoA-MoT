import collections

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell."""
    r, c = cell
    return {
        (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
        (r,     c - 1),             (r,     c + 1),
        (r + 1, c - 1), (r + 1, c), (r + 1, c + 1),
    }

def simulate_step(live_cells):
    """Simulates one step of Conway's Game of Life."""
    if not live_cells:
        return set()

    # Find all cells that need to be considered: live cells and their neighbors
    potential_cells = collections.defaultdict(int)
    for cell in live_cells:
        for neighbor in get_neighbors(cell):
            potential_cells[neighbor] += 1

    next_generation = set()
    for cell, count in potential_cells.items():
        # Rule for reproduction: a dead cell with exactly 3 live neighbors becomes alive
        if cell not in live_cells and count == 3:
            next_generation.add(cell)
        # Rule for survival: a live cell with 2 or 3 live neighbors survives
        elif cell in live_cells and (count == 2 or count == 3):
            next_generation.add(cell)
    
    return next_generation

def normalize_pattern(live_cells):
    """Normalizes a pattern's coordinates to detect cycles regardless of position."""
    if not live_cells:
        return frozenset()
    
    min_r = min(r for r, c in live_cells)
    min_c = min(c for r, c in live_cells)
    
    return frozenset((r - min_r, c - min_c) for r, c in live_cells)

def main():
    """
    Calculates how many of the 512 possible 3x3 starting grids in Conway's Game of Life
    eventually result in no living cells.
    """
    total_configs = 2**9
    dying_configs_count = 0
    max_generations = 200 # A safe upper limit to prevent any unforeseen infinite loops

    for i in range(total_configs):
        # Create the initial 3x3 grid from the integer i
        live_cells = set()
        for cell_index in range(9):
            if (i >> cell_index) & 1:
                row = cell_index // 3
                col = cell_index % 3
                live_cells.add((row, col))
        
        history = set()
        current_cells = live_cells
        
        for _ in range(max_generations):
            # Condition 1: The grid is empty
            if not current_cells:
                dying_configs_count += 1
                break

            # Condition 2: A cycle is detected
            normalized = normalize_pattern(current_cells)
            if normalized in history:
                break # This pattern will not die out
            
            history.add(normalized)
            
            # Evolve to the next generation
            current_cells = simulate_step(current_cells)

    print(f"Out of {total_configs} possible 3x3 configurations, {dying_configs_count} will eventually result in no living cells.")

if __name__ == "__main__":
    main()
<<<386>>>