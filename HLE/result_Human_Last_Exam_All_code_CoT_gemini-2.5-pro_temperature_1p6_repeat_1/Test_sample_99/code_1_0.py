import sys
from collections import Counter

def create_pn_pattern(n):
    """Creates the initial set of live cells for a Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, -i))
        live_cells.add((i, -i))
        live_cells.add((-i, i))
    return live_cells

def get_bounding_box_dims(live_cells):
    """Calculates the width and height of the bounding box of the pattern."""
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def next_generation(live_cells):
    """Computes the next generation of live cells based on Conway's rules."""
    neighbour_counts = Counter()
    for x, y in live_cells:
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j == 0:
                    continue
                neighbour_counts[(x + i, y + j)] += 1
    
    next_gen_cells = set()
    for cell, count in neighbour_counts.items():
        # A dead cell with 3 neighbours becomes live
        if cell not in live_cells and count == 3:
            next_gen_cells.add(cell)
        # A live cell with 2 or 3 neighbours survives
        elif cell in live_cells and (count == 2 or count == 3):
            next_gen_cells.add(cell)
            
    return next_gen_cells

def solve():
    """
    Finds the smallest value for Pn > 0 that causes the pattern to grow
    to at least twice its original size along any dimension.
    """
    n = 0
    MAX_GENERATIONS = 500  # A safe upper limit for generations to simulate

    while True:
        n += 1
        initial_cells = create_pn_pattern(n)
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        current_cells = initial_cells
        # Use frozenset for a hashable set to store in history
        history = {frozenset(current_cells)}
        
        grew_to_size = False
        for _ in range(MAX_GENERATIONS):
            current_cells = next_generation(current_cells)
            
            # Check for extinction (pattern dies out)
            if not current_cells:
                break
                
            # Check for a stable state or a repeating cycle
            frozen_current_cells = frozenset(current_cells)
            if frozen_current_cells in history:
                break
            history.add(frozen_current_cells)

            # Check for the growth condition
            width, height = get_bounding_box_dims(current_cells)
            if width >= target_dimension or height >= target_dimension:
                # The question asks for the smallest value for Pn.
                # Since we are iterating n from 1, the first n that
                # satisfies the condition is the answer.
                print(n)
                grew_to_size = True
                break
        
        if grew_to_size:
            break

solve()
<<<5>>>