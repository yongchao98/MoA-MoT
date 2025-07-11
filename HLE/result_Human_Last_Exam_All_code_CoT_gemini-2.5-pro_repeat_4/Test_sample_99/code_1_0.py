def solve_game_of_life_pn():
    """
    Finds the smallest value for Pn > 0 that causes the pattern to grow
    to at least twice its original size along any dimension in Conway's Game of Life.
    """

    def create_pn_pattern(n):
        """Creates the set of live cells for the Pn pattern."""
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))
        return live_cells

    def step(live_cells):
        """Computes the next generation of live cells based on Conway's rules."""
        neighbor_counts = {}
        # Count neighbors for all cells adjacent to any live cell.
        for r, c in live_cells:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    neighbor = (r + dr, c + dc)
                    neighbor_counts[neighbor] = neighbor_counts.get(neighbor, 0) + 1

        next_live_cells = set()
        # Apply the rules:
        # 1. A dead cell with 3 live neighbors becomes a live cell.
        # 2. A live cell with 2 or 3 live neighbors stays alive.
        for cell, count in neighbor_counts.items():
            if count == 3 or (count == 2 and cell in live_cells):
                next_live_cells.add(cell)
        return next_live_cells

    def get_dimensions(live_cells):
        """Calculates the bounding box width and height of the pattern."""
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
        live_cells = create_pn_pattern(n)
        
        # Calculate initial and target dimensions.
        original_dim = 2 * n + 1
        target_dim = 2 * original_dim

        # Use a history of states to detect if the pattern stabilizes or enters a loop.
        history = {frozenset(live_cells)}
        
        # Cap the simulation at 500 generations for each n.
        max_generations = 500 

        for generation in range(max_generations):
            live_cells = step(live_cells)
            
            # If the pattern dies out or stabilizes, it fails for this n.
            if not live_cells or frozenset(live_cells) in history:
                break
            history.add(frozenset(live_cells))

            # Check if the pattern has grown to the target size.
            width, height = get_dimensions(live_cells)
            if width >= target_dim or height >= target_dim:
                # Found the answer, print the explanation.
                print("Searching for the smallest n > 0...")
                print(f"For Pn with n = {n}:")
                print(f"Initial size = (2 * {n}) + 1 = {original_dim}")
                print(f"Target size >= 2 * {original_dim} = {target_dim}")
                print(f"After {generation + 1} generations, the pattern grew to a size of {width}x{height}.")
                if width >= target_dim:
                    print(f"The final width {width} is >= the target dimension {target_dim}.")
                if height >= target_dim:
                    print(f"The final height {height} is >= the target dimension {target_dim}.")
                print(f"Thus, the smallest n is {n}.")
                return n

# Run the search function.
result = solve_game_of_life_pn()
print(f"<<<{result}>>>")