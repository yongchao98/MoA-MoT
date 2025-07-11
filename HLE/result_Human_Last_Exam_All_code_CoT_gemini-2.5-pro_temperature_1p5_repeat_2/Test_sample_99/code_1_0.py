import sys

def simulate_game_of_life_for_pn():
    """
    Finds the smallest integer n > 0 for which the Pn pattern in Conway's Game of Life
    grows to at least twice its original size in any dimension.
    """

    def get_neighbors(cell):
        """Returns the 8 neighbors of a given cell."""
        r, c = cell
        return {
            (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
            (r,     c - 1),             (r,     c + 1),
            (r + 1, c - 1), (r + 1, c), (r + 1, c + 1),
        }

    def step(live_cells):
        """Calculates the next generation of live cells."""
        if not live_cells:
            return set()

        # Identify all cells that could possibly be live in the next generation
        # (live cells and their neighbors)
        candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
        
        next_live_cells = set()
        for cell in candidates:
            live_neighbors_count = len(get_neighbors(cell).intersection(live_cells))
            
            # Rule 1 & 2: A live cell with 2 or 3 live neighbors survives
            if cell in live_cells and live_neighbors_count in [2, 3]:
                next_live_cells.add(cell)
            # Rule 4: A dead cell with exactly 3 live neighbors becomes a live cell
            elif cell not in live_cells and live_neighbors_count == 3:
                next_live_cells.add(cell)
                
        return next_live_cells

    def create_pn(n):
        """Creates the initial set of live cells for the Pn pattern."""
        cells = {(0, 0)}
        for i in range(1, n + 1):
            cells.add((i, i))
            cells.add((-i, i))
            cells.add((i, -i))
            cells.add((-i, -i))
        return cells

    def get_dimensions(live_cells):
        """Calculates the width and height of the bounding box of live cells."""
        if not live_cells:
            return 0, 0
        min_r = min(cell[0] for cell in live_cells)
        max_r = max(cell[0] for cell in live_cells)
        min_c = min(cell[1] for cell in live_cells)
        max_c = max(cell[1] for cell in live_cells)
        height = max_r - min_r + 1
        width = max_c - min_c + 1
        return width, height

    n = 0
    # Set a reasonable limit for generations to avoid true infinite loops
    max_generations = 300 

    while True:
        n += 1
        live_cells = create_pn(n)
        
        initial_dimension = 2 * n + 1
        target_dimension = 2 * initial_dimension
        
        # History to detect stable/oscillating patterns
        history = [live_cells]

        for gen in range(1, max_generations + 1):
            live_cells = step(live_cells)
            
            # Check for termination conditions
            if not live_cells or live_cells in history:
                break
            
            # Keep history short to save memory
            history.append(live_cells)
            if len(history) > 20:
                history.pop(0)

            width, height = get_dimensions(live_cells)
            
            # Check if growth condition is met
            if width >= target_dimension or height >= target_dimension:
                final_dimension = max(width, height)
                print(f"The smallest integer n > 0 is {n}.")
                print(f"For n = {n}:")
                print(f"The initial dimension is 2 * {n} + 1 = {initial_dimension}.")
                print(f"The target dimension is 2 * {initial_dimension} = {target_dimension}.")
                print(f"After {gen} generations, the pattern's bounding box is {width}x{height}.")
                print(f"Since {final_dimension} >= {target_dimension}, the condition is met.")
                # Found the smallest n, so we can exit
                return n

# Run the simulation and print the result
final_n = simulate_game_of_life_for_pn()
# Final answer in the specified format is extracted from the simulation result
sys.stdout.write(f"\n<<<{final_n}>>>\n")
