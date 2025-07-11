import collections

def generate_pn_cells(n):
    """Generates the initial set of live cells for the Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_pattern_dimensions(live_cells):
    """Calculates the width and height of the bounding box of live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(x for x, y in live_cells)
    max_x = max(x for x, y in live_cells)
    min_y = min(y for x, y in live_cells)
    max_y = max(y for x, y in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def run_game_of_life_step(live_cells):
    """Performs one step (generation) of Conway's Game of Life."""
    if not live_cells:
        return set()

    # Use a defaultdict to count neighbors of all potentially active cells.
    # A potentially active cell is any neighbor of a currently live cell.
    neighbor_counts = collections.defaultdict(int)
    for x, y in live_cells:
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                neighbor_counts[(x + dx, y + dy)] += 1

    next_generation_cells = set()
    for cell, count in neighbor_counts.items():
        # Rule for birth: A dead cell with exactly 3 live neighbors becomes a live cell.
        if cell not in live_cells and count == 3:
            next_generation_cells.add(cell)
        # Rule for survival: A live cell with 2 or 3 live neighbors survives.
        elif cell in live_cells and (count == 2 or count == 3):
            next_generation_cells.add(cell)
            
    return next_generation_cells

def find_smallest_growing_pn():
    """
    Finds the smallest n > 0 for which the Pn pattern grows to at least
    twice its original size in any dimension.
    """
    n = 0
    # Set a high enough limit to allow for slow-growing patterns.
    max_generations_to_check = 500

    while True:
        n += 1
        
        initial_cells = generate_pn_cells(n)
        
        # The initial dimension (width and height) of Pn is 2n+1.
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        current_cells = initial_cells
        # Keep a history of cell states to detect cycles.
        history = {frozenset(current_cells)}

        for gen in range(1, max_generations_to_check + 1):
            current_cells = run_game_of_life_step(current_cells)

            # Condition 1: Pattern died out. It won't grow.
            if not current_cells:
                break

            # Condition 2: Pattern entered a repeating cycle. It won't grow further.
            frozen_state = frozenset(current_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)
            
            # Condition 3: Pattern grew to the target size. We found the solution.
            current_width, current_height = get_pattern_dimensions(current_cells)
            if current_width >= target_dimension or current_height >= target_dimension:
                achieved_dimension = max(current_width, current_height)
                
                print(f"Found the smallest value for n > 0 is {n}.")
                print(f"The initial dimension for P{n} is given by the equation: 2 * {n} + 1 = {initial_dimension}.")
                print(f"The target dimension is twice the initial, from the equation: 2 * {initial_dimension} = {target_dimension}.")
                print(f"At generation {gen}, the pattern reached a dimension of {achieved_dimension}, which meets the target.")
                print(f"<<<{n}>>>")
                return

        # If the loop completes without meeting the growth condition, continue to the next n.

if __name__ == '__main__':
    find_smallest_growing_pn()