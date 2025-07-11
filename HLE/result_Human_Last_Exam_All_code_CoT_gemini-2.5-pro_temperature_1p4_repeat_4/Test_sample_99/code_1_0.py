import sys

def get_next_generation(live_cells):
    """Calculates the next state of the game based on Conway's rules."""
    # A cell's neighbors are the 8 cells surrounding it.
    # We only need to consider cells that are currently alive or are neighbors of live cells.
    potential_cells = set()
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                potential_cells.add((x + dx, y + dy))

    next_live_cells = set()
    for x, y in potential_cells:
        # Count live neighbors
        live_neighbors = 0
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                if (x + dx, y + dy) in live_cells:
                    live_neighbors += 1
        
        # Apply rules:
        # 1. A live cell with 2 or 3 live neighbors survives.
        # 2. A dead cell with exactly 3 live neighbors becomes a live cell.
        is_alive = (x, y) in live_cells
        if is_alive and live_neighbors in [2, 3]:
            next_live_cells.add((x, y))
        elif not is_alive and live_neighbors == 3:
            next_live_cells.add((x, y))
            
    return next_live_cells

def get_dimensions(live_cells):
    """Calculates the width and height of the bounding box around live cells."""
    if not live_cells:
        return 0, 0
    
    min_x = min(cell[0] for cell in live_cells)
    max_x = max(cell[0] for cell in live_cells)
    min_y = min(cell[1] for cell in live_cells)
    max_y = max(cell[1] for cell in live_cells)
    
    return (max_x - min_x + 1), (max_y - min_y + 1)

def create_pn_pattern(n):
    """Creates the initial set of live cells for a Pn pattern."""
    cells = {(0, 0)}
    for i in range(1, n + 1):
        cells.add((i, i))
        cells.add((-i, i))
        cells.add((i, -i))
        cells.add((-i, -i))
    return cells

def find_smallest_growing_pn():
    """
    Iterates through n=1, 2, 3... to find the first Pn pattern that
    grows to at least twice its original dimension.
    """
    n = 0
    # Search up to a reasonable limit. The actual answer is small.
    while n < 15:
        n += 1
        print(f"--- Checking Pn for n = {n} ---")
        live_cells = create_pn_pattern(n)
        
        original_dim = 2 * n + 1
        target_dim = 2 * original_dim
        
        history = {frozenset(live_cells)}
        
        # Simulate for a limited number of generations to find a result
        for gen in range(1, 500):
            live_cells = get_next_generation(live_cells)
            
            # Stop if pattern dies or stabilizes/oscillates
            if not live_cells:
                print(f"Result for n={n}: Pattern died out.")
                break
            
            state = frozenset(live_cells)
            if state in history:
                print(f"Result for n={n}: Pattern stabilized or entered a cycle.")
                break
            history.add(state)
            
            width, height = get_dimensions(live_cells)
            current_max_dim = max(width, height)
            
            # Check if the growth condition is met
            if current_max_dim >= target_dim:
                print(f"\nSOLUTION FOUND: The smallest value for n is {n}.")
                print("\nDetails of the final state:")
                print(f"The P{n} pattern's initial dimension was {original_dim}.")
                print(f"The target dimension was 2 * {original_dim} = {target_dim}.")
                print(f"At generation {gen}, the pattern's bounding box became {width}x{height}.")
                
                print(f"\nFinal Equation Check:")
                print(f"max({width}, {height}) >= 2 * (2 * {n} + 1)")
                print(f"  {current_max_dim} >= 2 * ({original_dim})")
                print(f"  {current_max_dim} >= {target_dim}")
                print(f"  This is {current_max_dim >= target_dim}.")
                return n
                
    return None

if __name__ == '__main__':
    solution = find_smallest_growing_pn()
    if solution is None:
        print("\nNo solution was found within the search limits.")
