import collections

def get_next_generation(live_cells):
    """
    Calculates the next generation of live cells in Conway's Game of Life.
    
    Args:
        live_cells: A set of (row, col) tuples representing live cells.
        
    Returns:
        A new set of (row, col) tuples for the next generation.
    """
    if not live_cells:
        return set()

    # Count neighbors for all relevant cells (live cells and their neighbors)
    neighbor_counts = collections.defaultdict(int)
    for r, c in live_cells:
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                neighbor_counts[(r + dr, c + dc)] += 1

    next_live_cells = set()
    # Apply birth rule: a dead cell with 3 neighbors becomes live
    for cell, count in neighbor_counts.items():
        if cell not in live_cells and count == 3:
            next_live_cells.add(cell)
    
    # Apply survival rule: a live cell with 2 or 3 neighbors survives
    for cell in live_cells:
        if neighbor_counts.get(cell, 0) in [2, 3]:
            next_live_cells.add(cell)
            
    return next_live_cells

def get_bounding_box_size(live_cells):
    """
    Calculates the width and height of the bounding box around live cells.
    """
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def find_smallest_growing_pn():
    """
    Finds the smallest n > 0 for which the Pn pattern grows to at least
    twice its original size in any dimension.
    """
    n = 0
    max_generations_to_check = 500

    # Loop through n = 1, 2, 3, ...
    while True:
        n += 1
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        # Create the initial Pn pattern
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))

        history = {frozenset(live_cells)}
        
        # Simulate generations for the current n
        for generation in range(1, max_generations_to_check + 1):
            live_cells = get_next_generation(live_cells)
            
            # Condition 1: Pattern died out
            if not live_cells:
                break
            
            # Condition 2: Pattern stabilized or entered a cycle
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)
            
            # Condition 3: Pattern grew to the target size
            width, height = get_bounding_box_size(live_cells)
            if width >= target_dimension or height >= target_dimension:
                print(f"The smallest value for n > 0 is {n}.")
                print("\nHere is the breakdown:")
                print(f"For the pattern Pn where n = {n}:")
                print(f"The initial dimension is calculated as: 2 * n + 1 = 2 * {n} + 1 = {initial_dimension}")
                print(f"The target dimension is twice the initial: 2 * {initial_dimension} = {target_dimension}")
                print(f"At generation {generation}, the pattern grew to a size of {width}x{height}, which meets the target.")
                return

if __name__ == '__main__':
    find_smallest_growing_pn()