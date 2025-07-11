import sys

def solve_game_of_life():
    """
    Finds the smallest integer n > 0 for which the Pn pattern in Conway's Game of Life
    grows to at least twice its original size in any dimension.
    """

    def create_p_n(n):
        """Creates the Pn pattern as a set of (x, y) coordinates."""
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))
            live_cells.add((-i, i))
        return live_cells

    def get_bounding_box_dims(live_cells):
        """Calculates the dimensions of the bounding box of the live cells."""
        if not live_cells:
            return 0, 0, (0, 0, 0, 0)
        min_x = min(c[0] for c in live_cells)
        max_x = max(c[0] for c in live_cells)
        min_y = min(c[1] for c in live_cells)
        max_y = max(c[1] for c in live_cells)
        dim_x = max_x - min_x + 1
        dim_y = max_y - min_y + 1
        box = (min_x, max_x, min_y, max_y)
        return dim_x, dim_y, box

    def step(live_cells):
        """Computes the next generation of live cells based on the rules."""
        # Candidates for next generation are live cells and their neighbors
        candidates = set()
        for x, y in live_cells:
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    candidates.add((x + dx, y + dy))

        next_gen = set()
        for x, y in candidates:
            # Count live neighbors
            neighbors = 0
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    if (x + dx, y + dy) in live_cells:
                        neighbors += 1
            
            # Apply Game of Life rules
            is_alive = (x, y) in live_cells
            if is_alive and neighbors in [2, 3]:
                next_gen.add((x, y))
            elif not is_alive and neighbors == 3:
                next_gen.add((x, y))
        return next_gen

    n = 0
    while True:
        n += 1
        live_cells = create_p_n(n)
        
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        history = {frozenset(live_cells)}
        
        # Set a reasonable generation limit to prevent any unforeseen infinite runs.
        # Patterns that produce spaceships will meet the criteria relatively quickly.
        max_generations = 500 

        for generation in range(max_generations):
            live_cells = step(live_cells)
            
            if not live_cells:
                # Pattern died out
                break

            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                # Pattern stabilized or entered a repeating cycle
                break
            history.add(frozen_state)

            current_dim_x, current_dim_y, box = get_bounding_box_dims(live_cells)

            if current_dim_x >= target_dim or current_dim_y >= target_dim:
                # Found the solution
                print(f"The smallest value for Pn is n = {n}.")
                print(f"For P{n}, the initial dimension is 2 * {n} + 1 = {initial_dim}.")
                print(f"The target dimension to reach is 2 * {initial_dim} = {target_dim}.")
                print(f"At generation {generation + 1}, the pattern's bounding box is from ({box[0]}, {box[2]}) to ({box[1]}, {box[3]}).")
                print(f"The size is {current_dim_x}x{current_dim_y}, which meets the condition.")
                print(f"Final Answer: {n}")
                return n

# Run the simulation to find the answer
solve_game_of_life()
<<<6>>>