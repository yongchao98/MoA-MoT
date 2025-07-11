import sys
from collections import Counter

def create_pn_pattern(n):
    """Creates the initial set of live cells for a Pn pattern."""
    live_cells = set([(0, 0)])
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_next_generation(live_cells):
    """Calculates the next state of the grid based on Conway's Game of Life rules."""
    if not live_cells:
        return set()

    # Use a Counter to find the number of live neighbors for each cell.
    # We only need to check cells that are neighbors of current live cells.
    neighbor_counts = Counter()
    for x, y in live_cells:
        for i in range(x - 1, x + 2):
            for j in range(y - 1, y + 2):
                if (i, j) != (x, y):
                    neighbor_counts[(i, j)] += 1

    next_live_cells = set()
    for cell, count in neighbor_counts.items():
        # A dead cell with 3 neighbors becomes alive.
        if count == 3 and cell not in live_cells:
            next_live_cells.add(cell)
        # A live cell with 2 or 3 neighbors survives.
        elif count in (2, 3) and cell in live_cells:
            next_live_cells.add(cell)
            
    return next_live_cells

def get_bounding_box_dims(live_cells):
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
    Finds the smallest integer n > 0 for which the Pn pattern grows
    to at least twice its original size along any dimension.
    """
    MAX_GENERATIONS = 5000
    n = 0
    while True:
        n += 1
        
        cells = create_pn_pattern(n)
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        for _ in range(MAX_GENERATIONS):
            cells = get_next_generation(cells)
            
            # If the pattern dies out, it can't grow. Move to the next n.
            if not cells:
                break
                
            width, height = get_bounding_box_dims(cells)
            
            # Check if the growth condition is met.
            if width >= target_dimension or height >= target_dimension:
                print(n)
                return

if __name__ == '__main__':
    solve()