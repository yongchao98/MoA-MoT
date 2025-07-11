import sys

def create_pn_pattern(n):
    """Creates the initial Pn pattern as a set of live cell coordinates."""
    cells = {(0, 0)}
    for i in range(1, n + 1):
        cells.add((i, i))
        cells.add((-i, -i))
        cells.add((i, -i))
        cells.add((-i, i))
    return cells

def get_neighbors(cell):
    """Returns a set of 8 neighbor coordinates for a given cell."""
    x, y = cell
    return {
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    }

def calculate_next_generation(live_cells):
    """
    Calculates the next state of the grid based on Conway's rules.
    - A live cell with 2 or 3 live neighbors survives.
    - A dead cell with exactly 3 live neighbors becomes a live cell.
    """
    # Consider all live cells and their neighbors as potentially live in the next generation.
    potential_cells = live_cells.union(set(neighbor for cell in live_cells for neighbor in get_neighbors(cell)))
    next_gen_cells = set()
    
    for cell in potential_cells:
        live_neighbor_count = len(get_neighbors(cell) & live_cells)
        if live_neighbor_count == 3 or (live_neighbor_count == 2 and cell in live_cells):
            next_gen_cells.add(cell)
            
    return next_gen_cells

def get_pattern_dimensions(live_cells):
    """Calculates the width and height of the pattern's bounding box."""
    if not live_cells:
        return 0, 0
    min_x = min(cell[0] for cell in live_cells)
    max_x = max(cell[0] for cell in live_cells)
    min_y = min(cell[1] for cell in live_cells)
    max_y = max(cell[1] for cell in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def solve():
    """
    Finds the smallest n > 0 for which Pn grows to at least twice its initial size.
    """
    n = 0
    # A sufficiently large number of generations to test for growth or stability.
    MAX_GENERATIONS = 500

    while True:
        n += 1
        live_cells = create_pn_pattern(n)
        
        # Define initial and target sizes
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        # Keep a history of states to detect stable cycles
        history = {frozenset(live_cells)}
        
        grew_large_enough = False
        
        for generation in range(MAX_GENERATIONS):
            live_cells = calculate_next_generation(live_cells)
            
            if not live_cells:  # Extinction
                break
                
            current_state = frozenset(live_cells)
            if current_state in history:  # Stable or oscillating pattern
                break
            history.add(current_state)
            
            width, height = get_pattern_dimensions(live_cells)
            if width >= target_dimension or height >= target_dimension:
                grew_large_enough = True
                break
                
        if grew_large_enough:
            # We found the smallest n that satisfies the condition.
            # The prompt requires printing the final number.
            print(n)
            return

if __name__ == '__main__':
    solve()