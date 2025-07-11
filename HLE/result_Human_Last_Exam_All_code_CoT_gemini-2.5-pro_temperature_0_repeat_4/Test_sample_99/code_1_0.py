def solve_conway_pn():
    """
    In Conway's Game of Life, Pn denotes an X-shaped starting pattern such that
    a central live cell has n live squares radiating in a diagonal line from
    each of its corners. This script finds the smallest value for n > 0 that
    causes the pattern to grow to at least twice its original size along any dimension.
    """

    def create_pattern_pn(n):
        """Creates the set of live cells for a Pn pattern."""
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))
        return live_cells

    def evolve(live_cells):
        """Evolves the set of live cells by one generation."""
        # Find all cells to consider (live cells + their neighbors)
        candidates = set()
        for (x, y) in live_cells:
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    candidates.add((x + dx, y + dy))

        next_gen_cells = set()
        for (x, y) in candidates:
            # Count live neighbors
            neighbors = 0
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    if (x + dx, y + dy) in live_cells:
                        neighbors += 1

            # Apply Conway's Game of Life rules
            is_alive = (x, y) in live_cells
            if is_alive and (neighbors == 2 or neighbors == 3):
                next_gen_cells.add((x, y))
            elif not is_alive and neighbors == 3:
                next_gen_cells.add((x, y))

        return next_gen_cells

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

    # Set reasonable limits for the search
    MAX_N = 20
    MAX_GENERATIONS = 1000

    # Iterate through n = 1, 2, 3, ...
    for n in range(1, MAX_N + 1):
        live_cells = create_pattern_pn(n)
        
        # The initial pattern is symmetric, so width and height are the same.
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        history = set()
        for generation in range(MAX_GENERATIONS):
            # Detect stable or oscillating patterns to avoid infinite loops
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break  # Pattern is stable/oscillating, try next n
            history.add(frozen_state)

            live_cells = evolve(live_cells)

            if not live_cells:
                break # Pattern died out, try next n

            # Check if the growth condition is met
            current_width, current_height = get_dimensions(live_cells)
            if current_width >= target_dim or current_height >= target_dim:
                # Found the solution. Print the value of n.
                # The "equation" is that the final dimension is at least twice the initial.
                # The numbers in this relationship are n, the final dimension, and the initial dimension.
                print(n)
                return

if __name__ == '__main__':
    solve_conway_pn()