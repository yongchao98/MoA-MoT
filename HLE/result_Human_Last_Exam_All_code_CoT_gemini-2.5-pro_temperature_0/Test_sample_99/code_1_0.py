import collections

def solve_game_of_life_problem():
    """
    Finds the smallest integer n > 0 for which the Pn pattern in Conway's Game of Life
    grows to at least twice its original size in any dimension.
    """

    def get_next_generation(live_cells):
        """Calculates the next state of live cells in the Game of Life."""
        if not live_cells:
            return set()

        # Use a Counter to efficiently find all cells with live neighbors
        neighbor_counts = collections.Counter()
        for x, y in live_cells:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if i == 0 and j == 0:
                        continue
                    neighbor_counts[(x + i, y + j)] += 1

        next_live_cells = set()
        # Apply the rules of the game
        for cell, count in neighbor_counts.items():
            # A dead cell with 3 neighbors becomes live
            if count == 3 and cell not in live_cells:
                next_live_cells.add(cell)
            # A live cell with 2 or 3 neighbors survives
            elif count in (2, 3) and cell in live_cells:
                next_live_cells.add(cell)

        return next_live_cells

    def get_bounding_box_dims(live_cells):
        """Calculates the width and height of the bounding box for live cells."""
        if not live_cells:
            return 0, 0

        min_x = min(c[0] for c in live_cells)
        max_x = max(c[0] for c in live_cells)
        min_y = min(c[1] for c in live_cells)
        max_y = max(c[1] for c in live_cells)

        width = max_x - min_x + 1
        height = max_y - min_y + 1
        return width, height

    n = 0
    while True:
        n += 1
        
        # 1. Create the initial Pn pattern
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))

        # 2. Calculate initial and target dimensions
        initial_width, initial_height = get_bounding_box_dims(live_cells)
        # Since the pattern is symmetric, width and height are the same
        target_dimension = 2 * initial_width

        # 3. Simulate generations
        history = {frozenset(live_cells)}
        max_generations = 200  # A safe limit to detect growth

        current_cells = live_cells
        for gen in range(1, max_generations + 1):
            current_cells = get_next_generation(current_cells)

            # Check for termination conditions (pattern died or stabilized)
            if not current_cells:
                break
            
            frozen_state = frozenset(current_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            # 4. Check if the growth condition is met
            current_width, current_height = get_bounding_box_dims(current_cells)
            if current_width >= target_dimension or current_height >= target_dimension:
                print(f"Solution found for Pn where n = {n}.")
                print(f"The initial pattern P{n} has {4*n + 1} cells.")
                print(f"The initial dimension is {initial_width}x{initial_height}.")
                print(f"The target dimension was to reach at least {target_dimension}.")
                print(f"After {gen} generations, the pattern grew to {current_width}x{current_height}.")
                
                # Output the final equation check as requested
                if current_width >= target_dimension:
                    print(f"The width check is: {current_width} >= 2 * {initial_width}, which is True.")
                if current_height >= target_dimension:
                     print(f"The height check is: {current_height} >= 2 * {initial_height}, which is True.")
                
                return n

# Run the simulation and print the final answer
smallest_n = solve_game_of_life_problem()
print(f"\nThe smallest value for n > 0 is {smallest_n}.")
print(f"<<<{smallest_n}>>>")
