import collections

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
    """Calculates the next state of the grid based on Conway's rules."""
    # Using a Counter is an efficient way to find all cells that are
    # neighbors of live cells and count their number of live neighbors.
    neighbor_counts = collections.Counter()
    for r, c in live_cells:
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                neighbor = (r + dr, c + dc)
                neighbor_counts[neighbor] += 1

    next_live_cells = set()
    # Apply the rules of life.
    for cell, count in neighbor_counts.items():
        # A dead cell with exactly 3 live neighbors becomes a live cell.
        if count == 3 and cell not in live_cells:
            next_live_cells.add(cell)
        # A live cell with 2 or 3 live neighbors survives.
        elif count in (2, 3) and cell in live_cells:
            next_live_cells.add(cell)
            
    return next_live_cells

def get_bounding_box_size(live_cells):
    """Calculates the width and height of the pattern's bounding box."""
    if not live_cells:
        return 0, 0
    min_r = min(r for r, c in live_cells)
    max_r = max(r for r, c in live_cells)
    min_c = min(c for r, c in live_cells)
    max_c = max(c for r, c in live_cells)
    
    height = max_r - min_r + 1
    width = max_c - min_c + 1
    return width, height

def solve():
    """
    Finds the smallest n > 0 for Pn that grows to at least
    twice its original size in any dimension.
    """
    n = 0
    max_generations_per_n = 200 # A safe limit to prevent infinite loops

    while True:
        n += 1
        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        live_cells = create_pn_pattern(n)
        
        # Keep a history to detect stable patterns or cycles.
        history = {frozenset(live_cells)}

        for gen in range(max_generations_per_n):
            live_cells = get_next_generation(live_cells)

            # Condition 1: Pattern dies out.
            if not live_cells:
                break

            # Condition 2: Pattern stabilizes or enters a cycle.
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            # Condition 3: Check if the pattern has grown enough.
            width, height = get_bounding_box_size(live_cells)
            if width >= target_size:
                print(f"{width} >= 2 * (2 * {n} + 1)")
                return n
            if height >= target_size:
                print(f"{height} >= 2 * (2 * {n} + 1)")
                return n
    return -1 # Should not be reached if a solution exists

if __name__ == '__main__':
    result = solve()
    if result != -1:
        print(f'<<<{result}>>>')
