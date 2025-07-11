import sys
from collections import Counter

def create_pn_pattern(n):
    """Generates the Pn pattern, an X-shape with arms of length n."""
    cells = set([(0, 0)])
    for i in range(1, n + 1):
        cells.add((i, i))
        cells.add((-i, i))
        cells.add((i, -i))
        cells.add((-i, -i))
    return cells

def step(live_cells):
    """
    Calculates the next generation of live cells based on Conway's rules.
    1. A live cell with 2 or 3 live neighbours survives.
    2. A dead cell with 3 live neighbours becomes a live cell.
    """
    if not live_cells:
        return set()
    
    # Count live neighbours for all relevant cells (live cells and their neighbours)
    neighbour_counts = Counter()
    for x, y in live_cells:
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                neighbour_counts[(x + dx, y + dy)] += 1
    
    next_gen_cells = set()
    for cell, count in neighbour_counts.items():
        if count == 3 or (cell in live_cells and count == 2):
            next_gen_cells.add(cell)
            
    return next_gen_cells

def get_dimensions(live_cells):
    """Calculates the width and height of the pattern's bounding box."""
    if not live_cells:
        return 0, 0
    
    min_x = min(x for x, y in live_cells)
    max_x = max(x for x, y in live_cells)
    min_y = min(y for x, y in live_cells)
    max_y = max(y for x, y in live_cells)
    
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def solve_game_of_life_problem():
    """
    Finds the smallest integer n > 0 for which the Pn pattern grows
    to at least twice its original size in any dimension.
    """
    MAX_GENERATIONS = 300 
    
    n = 1
    while True:
        live_cells = create_pn_pattern(n)
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        history = set()

        for generation in range(MAX_GENERATIONS):
            width, height = get_dimensions(live_cells)
            
            if width >= target_dim or height >= target_dim:
                # Found the solution
                print(f"The smallest value for n is {n}.")
                print(f"Initial dimension for P{n}: {initial_dim}x{initial_dim}")
                print(f"Target dimension: {target_dim}")
                print(f"Dimension after {generation} generations: {width}x{height}")
                return n

            # Use frozenset for efficient cycle detection
            pattern_state = frozenset(live_cells)
            if pattern_state in history:
                # Cycle detected, pattern will not grow further for this n.
                break
            history.add(pattern_state)

            live_cells = step(live_cells)
            
            if not live_cells:
                # Pattern died out.
                break
        
        # If loop finishes without finding a solution for this n, try the next one.
        n += 1

if __name__ == '__main__':
    solve_game_of_life_problem()
