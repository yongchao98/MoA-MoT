import collections

def generate_pn(n):
    """Generates the set of live cells for the Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_bounding_box_dims(live_cells):
    """Calculates the width and height of the bounding box for the live cells."""
    if not live_cells:
        return 0, 0
    
    min_x = min(cell[0] for cell in live_cells)
    max_x = max(cell[0] for cell in live_cells)
    min_y = min(cell[1] for cell in live_cells)
    max_y = max(cell[1] for cell in live_cells)
    
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def get_next_generation(live_cells):
    """Calculates the next state of the cells based on Conway's rules."""
    if not live_cells:
        return set()
    
    # Count live neighbors for all relevant cells
    neighbor_counts = collections.Counter()
    for cell in live_cells:
        x, y = cell
        for i in range(x - 1, x + 2):
            for j in range(y - 1, y + 2):
                if (i, j) != cell:
                    neighbor_counts[(i, j)] += 1
            
    next_live_cells = set()
    # Apply the rules of life
    for cell, count in neighbor_counts.items():
        # A dead cell with 3 neighbors becomes live
        if cell not in live_cells and count == 3:
            next_live_cells.add(cell)
        # A live cell with 2 or 3 neighbors survives
        elif cell in live_cells and (count == 2 or count == 3):
            next_live_cells.add(cell)
            
    return next_live_cells

def solve_game_of_life_problem():
    """
    Finds the smallest n > 0 for which Pn grows to at least twice its
    original size in any dimension.
    """
    n = 0
    max_n_to_check = 10  # A reasonable upper limit for the search
    max_generations = 500 # A safety limit for each simulation

    while n < max_n_to_check:
        n += 1
        print(f"--- Checking Pn for n = {n} ---")
        
        live_cells = generate_pn(n)
        initial_width, initial_height = get_bounding_box_dims(live_cells)
        target_dimension = 2 * initial_width  # Since width and height are the same

        # Use frozenset for history check as sets are not hashable
        history = {frozenset(live_cells)}

        for gen in range(max_generations):
            live_cells = get_next_generation(live_cells)
            
            # Condition 1: Pattern died out
            if not live_cells:
                print(f"Result: Pattern died out at generation {gen + 1}.")
                break

            # Condition 2: Pattern stabilized or is oscillating
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                print(f"Result: Pattern stabilized or entered a loop at generation {gen + 1}.")
                break
            history.add(frozen_state)

            # Condition 3: Check for required growth
            current_width, current_height = get_bounding_box_dims(live_cells)
            if current_width >= target_dimension or current_height >= target_dimension:
                print(f"Success! Found the solution for Pn where n = {n}.")
                print(f"Initial size was {initial_width}x{initial_height}.")
                print(f"The target dimension was >= {target_dimension}.")
                print(f"At generation {gen + 1}, the pattern reached a size of {current_width}x{current_height}.")
                print(f"\nThe smallest value for Pn > 0 that causes the pattern to grow to at least twice its original size is {n}.")
                return n
        
        if len(live_cells) > 0:
             print(f"Result: Condition not met within {max_generations} generations.")

    print("Could not find a solution within the checked range of n.")
    return None

# Run the solver
solve_game_of_life_problem()