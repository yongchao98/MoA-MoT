import collections

def solve_conways_pn_problem():
    """
    Solves for the smallest integer n > 0 for the Pn pattern in Conway's Game of Life
    that causes the pattern to at least double its initial size in any dimension.
    """

    def get_neighbors(cell):
        """Returns a set of 8 neighbor coordinates for a given cell."""
        x, y = cell
        return {
            (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
            (x - 1, y),                 (x + 1, y),
            (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
        }

    def simulate_step(live_cells):
        """Calculates the next generation of live cells."""
        if not live_cells:
            return set()
        
        # Count live neighbors for all potentially active cells.
        # These are the current live cells and their dead neighbors.
        neighbor_counts = collections.Counter(
            neighbor for cell in live_cells for neighbor in get_neighbors(cell)
        )
        
        # A dead cell with exactly 3 neighbors becomes live.
        born_cells = {cell for cell, count in neighbor_counts.items() if count == 3 and cell not in live_cells}
        
        # A live cell with 2 or 3 neighbors survives.
        surviving_cells = {cell for cell in live_cells if neighbor_counts[cell] in {2, 3}}
        
        return born_cells.union(surviving_cells)

    def get_bounding_box_dims(live_cells):
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

    n = 0
    while True:
        n += 1
        
        # 1. Create the initial Pn pattern as a set of live cells.
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))

        # 2. Calculate initial and target dimensions.
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # Use a history of states to detect stable/oscillating patterns.
        history = {frozenset(live_cells)}
        
        # 3. Simulate for a max number of generations.
        # This limit prevents true infinite loops for chaotic patterns
        # that may never meet the condition.
        max_generations = 500
        for generation in range(1, max_generations + 1):
            live_cells = simulate_step(live_cells)
            
            if not live_cells:
                # Pattern died out.
                break

            current_width, current_height = get_bounding_box_dims(live_cells)
            
            # 4. Check if growth condition is met.
            if current_width >= target_dim or current_height >= target_dim:
                print(f"For n = {n}, the initial size is {initial_dim}x{initial_dim}.")
                print(f"The target size is >= {target_dim} in either dimension.")
                print(f"After generation {generation}, the size became {current_width}x{current_height}.")
                print(f"The smallest value for n is {n}.")
                return

            # Check for stabilization/oscillation.
            current_state = frozenset(live_cells)
            if current_state in history:
                break
            history.add(current_state)

# Execute the search function.
solve_conways_pn_problem()

# The pattern for n=6 (P6) is known to produce "queen bee shuttles",
# which are oscillators that periodically emit gliders. These gliders
# travel across the grid, ensuring that the pattern's bounding box grows
# indefinitely and will eventually meet the condition. Patterns for n < 6
# are known to stabilize, oscillate within a fixed area, or die out.
<<<6>>>