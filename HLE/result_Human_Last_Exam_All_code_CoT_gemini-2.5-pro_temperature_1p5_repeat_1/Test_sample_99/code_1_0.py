import sys

def create_pn(n):
    """Creates the initial set of live cells for a Pn pattern."""
    cells = set()
    cells.add((0, 0))
    for i in range(1, n + 1):
        cells.add((i, i))
        cells.add((i, -i))
        cells.add((-i, i))
        cells.add((-i, -i))
    return cells

def get_dimensions(live_cells):
    """Calculates the width and height of the bounding box of live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def step(live_cells):
    """Performs one generation step of Conway's Game of Life."""
    # A candidate for changing state is any live cell or any of its neighbors.
    candidates = set(live_cells)
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx != 0 or dy != 0:
                    candidates.add((x + dx, y + dy))
    
    next_live_cells = set()
    for x, y in candidates:
        # Count live neighbors
        count = 0
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if (dx != 0 or dy != 0) and (x + dx, y + dy) in live_cells:
                    count += 1
        
        # Apply the rules of life
        is_live = (x, y) in live_cells
        if is_live and count in [2, 3]:
            next_live_cells.add((x, y))
        elif not is_live and count == 3:
            next_live_cells.add((x, y))
            
    return next_live_cells

def solve():
    """
    Finds the smallest n > 0 for which Pn grows to at least twice its size.
    """
    n = 1
    # Set a practical limit for n and generations to prevent infinite loops.
    MAX_N = 20
    MAX_GENERATIONS = 500

    while n <= MAX_N:
        print(f"Testing Pn for n={n}...")
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim
        
        cells = create_pn(n)
        history = {frozenset(cells)}
        
        found_solution = False
        for gen in range(1, MAX_GENERATIONS + 1):
            cells = step(cells)
            
            # Condition 1: Pattern died
            if not cells:
                print(f"n={n}: Died out after {gen} generations.")
                break
            
            # Condition 2: Pattern is stable or an oscillator
            frozen_cells = frozenset(cells)
            if frozen_cells in history:
                print(f"n={n}: Stabilized after {gen} generations.")
                break
            history.add(frozen_cells)
            
            # Condition 3: Pattern has grown to the target size
            width, height = get_dimensions(cells)
            if width >= target_dim or height >= target_dim:
                print("\n" + "="*40)
                print("SOLUTION FOUND!")
                print("="*40)
                print(f"The smallest value for Pn is when n = {n}.")
                print("\nDetails of the solution:")
                print(f"Value of n: {n}")
                print(f"Initial Dimension: {initial_dim}")
                print(f"Target Dimension (>= 2 * Initial): {target_dim}")
                print(f"Generations to Grow: {gen}")
                print(f"Final Dimension: {width}x{height}")
                found_solution = True
                break
        
        if found_solution:
            return
        
        n += 1

    print("\nCould not find a solution within the defined limits.")

if __name__ == '__main__':
    solve()