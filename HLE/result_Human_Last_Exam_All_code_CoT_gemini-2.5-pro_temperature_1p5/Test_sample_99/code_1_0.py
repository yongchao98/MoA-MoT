import collections

def solve_game_of_life_problem():
    """
    Finds the smallest integer n > 0 for which the Pn pattern in Conway's Game of Life
    grows to at least twice its original size in any dimension.
    """

    def create_pn(n):
        """Creates the initial set of live cells for a Pn pattern."""
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))
        return live_cells

    def get_dimensions(live_cells):
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

    def step(live_cells):
        """Calculates the next generation of live cells."""
        neighbor_counts = collections.Counter()
        for x, y in live_cells:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if i == 0 and j == 0:
                        continue
                    neighbor_counts[(x + i, y + j)] += 1

        next_gen_live_cells = set()
        # Check all cells that are live or are neighbors of live cells
        all_potential_cells = set(live_cells).union(set(neighbor_counts.keys()))

        for cell in all_potential_cells:
            count = neighbor_counts.get(cell, 0)
            is_live = cell in live_cells
            
            # Rule 4: Reproduction
            if not is_live and count == 3:
                next_gen_live_cells.add(cell)
            # Rule 2 & 3: Survival
            elif is_live and (count == 2 or count == 3):
                next_gen_live_cells.add(cell)
            # Rule 1 (underpopulation) and Rule 4 (overpopulation) are handled by omission

        return next_gen_live_cells

    max_generations = 300
    n = 0
    while True:
        n += 1
        print(f"--- Testing n = {n} ---")
        
        initial_cells = create_pn(n)
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        current_cells = initial_cells
        history = {frozenset(current_cells)}
        
        found = False
        for gen in range(max_generations):
            current_cells = step(current_cells)
            width, height = get_dimensions(current_cells)

            # Check for growth condition
            if width >= target_dim or height >= target_dim:
                print(f"Success for n={n} at generation {gen + 1}.")
                print(f"Initial dimension was: 2 * {n} + 1 = {initial_dim}")
                print(f"Target dimension was: 2 * {initial_dim} = {target_dim}")
                final_dim = max(width, height)
                print(f"The pattern's dimension reached {final_dim}.")
                print(f"Final check: {final_dim} >= {target_dim}")
                print("\nSmallest value for n found.")
                print(f"Final Answer: {n}")
                return n

            # Check for extinction or stability/oscillation
            if not current_cells:
                print(f"Pattern for n={n} died out.")
                break
            
            frozen_cells = frozenset(current_cells)
            if frozen_cells in history:
                print(f"Pattern for n={n} stabilized or entered a cycle.")
                break
            history.add(frozen_cells)
        else:
             print(f"Pattern for n={n} did not reach target size in {max_generations} generations.")

# Execute the solver
solve_game_of_life_problem()