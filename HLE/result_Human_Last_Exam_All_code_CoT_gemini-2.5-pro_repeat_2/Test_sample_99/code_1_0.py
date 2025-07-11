import sys

def create_pn(n):
    """Creates the initial set of live cells for a Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_bounding_box_dims(live_cells):
    """Calculates the width and height of the bounding box of live cells."""
    if not live_cells:
        return 0, 0
    
    # sys.maxsize is used for finding the minimum coordinate value effectively.
    min_x = sys.maxsize
    max_x = -sys.maxsize
    min_y = sys.maxsize
    max_y = -sys.maxsize
    
    for x, y in live_cells:
        if x < min_x: min_x = x
        if x > max_x: max_x = x
        if y < min_y: min_y = y
        if y > max_y: max_y = y
        
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def step(live_cells):
    """Calculates the next generation of live cells based on Conway's rules."""
    if not live_cells:
        return set()

    neighbor_counts = {}
    for x, y in live_cells:
        # Count neighbors for all surrounding cells, including the cell itself
        for i in range(x - 1, x + 2):
            for j in range(y - 1, y + 2):
                if i == x and j == y:
                    continue
                neighbor_counts[(i, j)] = neighbor_counts.get((i, j), 0) + 1
    
    next_live_cells = set()
    for cell, count in neighbor_counts.items():
        # A dead cell with 3 neighbors becomes live
        if cell not in live_cells and count == 3:
            next_live_cells.add(cell)
        # A live cell with 2 or 3 neighbors survives
        elif cell in live_cells and (count == 2 or count == 3):
            next_live_cells.add(cell)
            
    return next_live_cells

def find_smallest_growing_pn():
    """
    Iterates through n=1, 2, 3... to find the smallest Pn that meets the growth criteria.
    """
    n = 0
    # Set a reasonable limit to avoid true infinite loops for complex patterns.
    max_generations = 500

    while True:
        n += 1
        live_cells = create_pn(n)
        
        # Pn patterns are symmetrical, so initial width and height are the same.
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # Store a history of states to detect stable or oscillating patterns.
        history = {frozenset(live_cells)}

        for generation in range(max_generations):
            live_cells = step(live_cells)
            
            # Stop if the pattern dies out
            if not live_cells:
                break
            
            # Stop if the pattern repeats (stabilizes or oscillates)
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            current_width, current_height = get_bounding_box_dims(live_cells)
            
            # Check if the growth condition is met
            if current_width >= target_dim or current_height >= target_dim:
                print(f"Solution found for Pn where n = {n}.")
                print(f"The initial pattern P{n} has dimensions {initial_dim}x{initial_dim}.")
                print(f"The growth condition is for either dimension to reach at least 2 * {initial_dim} = {target_dim}.")
                print(f"After {generation + 1} generations, the pattern's bounding box is {current_width}x{current_height}.")
                
                print("\nFinal equation check:")
                print(f"Width: {current_width} >= {target_dim} is {current_width >= target_dim}")
                print(f"Height: {current_height} >= {target_dim} is {current_height >= target_dim}")
                
                print(f"\nThe smallest value for n is {n}.")
                print(f"<<<{n}>>>")
                return

if __name__ == "__main__":
    find_smallest_growing_pn()