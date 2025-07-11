import sys

def get_neighbors(cell):
    """Returns the 8 neighbors of a given cell."""
    x, y = cell
    return [
        (x-1, y-1), (x, y-1), (x+1, y-1),
        (x-1, y),             (x+1, y),
        (x-1, y+1), (x, y+1), (x+1, y+1),
    ]

def step(live_cells):
    """Calculates the next generation of live cells."""
    next_live_cells = set()
    # A candidate for life is any current live cell or any of its neighbors.
    candidates = live_cells.union(set(n for cell in live_cells for n in get_neighbors(cell)))

    for cell in candidates:
        count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
        is_alive = cell in live_cells

        # Apply Conway's Game of Life rules
        if is_alive and count in (2, 3):
            next_live_cells.add(cell)
        elif not is_alive and count == 3:
            next_live_cells.add(cell)

    return next_live_cells

def create_pn(n):
    """Creates the initial set of live cells for a Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_dimensions(live_cells):
    """Calculates the width and height of the pattern's bounding box."""
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
    """Iterates through n to find the smallest value that meets the growth criteria."""
    n = 0
    # Set a reasonable limit on generations to check for growth, to avoid infinite loops.
    # For a pattern that generates gliders, this needs to be large enough for them
    # to travel sufficiently far.
    MAX_GENERATIONS = 500

    while True:
        n += 1
        live_cells = create_pn(n)

        initial_width, initial_height = get_dimensions(live_cells)
        target_size = 2 * initial_width # width and height are the same initially

        # Use a history to detect static patterns or oscillators
        history = {frozenset(live_cells)}
        
        current_cells = live_cells
        for gen in range(MAX_GENERATIONS):
            current_cells = step(current_cells)

            # Condition 1: Pattern died out
            if not current_cells:
                break
            
            # Condition 2: Pattern is static or oscillating
            frozen_set = frozenset(current_cells)
            if frozen_set in history:
                break
            history.add(frozen_set)
            
            # Condition 3: Pattern grew to at least twice its original size
            current_width, current_height = get_dimensions(current_cells)
            if current_width >= target_size or current_height >= target_size:
                # Found the solution, print the details
                print(f"Found the smallest value n = {n}.")
                print("\n--- Details ---")
                print("Equation for initial size: 2 * n + 1")
                print(f"Initial size for n={n}: 2 * {n} + 1 = {initial_width}")
                print("\nEquation for target size: 2 * initial_size")
                print(f"Target size: 2 * {initial_width} = {target_size}")

                print(f"\nAt generation {gen + 1}, the pattern grew to a size of {current_width}x{current_height}.")
                
                # Show which dimension met the condition
                if current_width >= target_size:
                    print(f"The width {current_width} is greater than or equal to the target size {target_size}.")
                if current_height >= target_size:
                     print(f"The height {current_height} is greater than or equal to the target size {target_size}.")

                return n # Return the answer

# Run the simulation and find the answer
result = find_smallest_growing_pn()
print(f"\n<<<__{result}__>>>")
