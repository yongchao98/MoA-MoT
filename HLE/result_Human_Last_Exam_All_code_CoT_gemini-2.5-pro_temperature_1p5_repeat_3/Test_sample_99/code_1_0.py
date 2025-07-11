import sys

def solve_game_of_life():
    """
    Finds the smallest integer n > 0 for pattern Pn that grows to at least 
    twice its original size in Conway's Game of Life.
    """
    n = 0
    found = False
    
    # The main loop iterates through n=1, 2, 3, ... until a solution is found.
    while not found:
        n += 1
        
        # --- Simulation setup for Pn ---
        # Create the initial Pn pattern as a set of (x, y) coordinates.
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))
            live_cells.add((-i, i))

        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        # History to detect stable patterns or non-moving oscillators.
        history = set()
        
        # We need a high enough generation limit as some patterns are "Methuselahs"
        # that evolve for a long time before stabilizing.
        max_generations = 2000 

        for g in range(max_generations):
            if not live_cells:
                # The pattern died out, so it cannot grow. Move to the next n.
                break

            # To detect still-lifes or simple oscillators, we store a snapshot
            # of the set of live cells. If a state repeats, we are in a loop.
            pattern_snapshot = frozenset(live_cells)
            if pattern_snapshot in history:
                break # Cycle detected. Move to the next n.
            history.add(pattern_snapshot)

            # --- Check the growth condition ---
            # Calculate the dimensions of the current bounding box.
            min_x = min(c[0] for c in live_cells)
            max_x = max(c[0] for c in live_cells)
            min_y = min(c[1] for c in live_cells)
            max_y = max(c[1] for c in live_cells)
            current_width = max_x - min_x + 1
            current_height = max_y - min_y + 1

            if current_width >= target_size or current_height >= target_size:
                # Success! We found the smallest n. Print the results.
                print(f"The smallest value for Pn > 0 is n = {n}.")
                print(f"The initial pattern P{n} has a bounding box of {initial_size}x{initial_size}.")
                print(f"The target dimension is at least twice the initial, which is >= 2 * {initial_size} = {target_size}.")
                print(f"After {g} generations, the pattern grows to a size of {current_width}x{current_height}.")
                
                # Output the final equation with the discovered numbers.
                final_dim = max(current_width, current_height)
                print(f"\nThe condition is met, as shown by the equation:")
                print(f"max(final_width, final_height) >= 2 * (2 * n + 1)")
                print(f"max({current_width}, {current_height}) >= 2 * (2 * {n} + 1)")
                print(f"{final_dim} >= {target_size}")
                
                # Set 'found' flag to True to exit the main while-loop.
                found = True
                break

            # --- Calculate the next generation ---
            next_live_cells = set()
            
            # To optimize, we only check cells that could possibly change state:
            # existing live cells and their immediate neighbors.
            candidates = set()
            for x, y in live_cells:
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        candidates.add((x + dx, y + dy))

            for x, y in candidates:
                # Count live neighbors for the candidate cell.
                neighbors = 0
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if dx == 0 and dy == 0:
                            continue
                        if (x + dx, y + dy) in live_cells:
                            neighbors += 1
                
                # Apply the standard rules of Conway's Game of Life.
                is_alive = (x, y) in live_cells
                if is_alive and (neighbors == 2 or neighbors == 3):
                    next_live_cells.add((x, y)) # Survival
                elif not is_alive and neighbors == 3:
                    next_live_cells.add((x, y)) # Birth
            
            live_cells = next_live_cells
        
    return n

# Execute the function and capture the final answer.
final_answer = solve_game_of_life()
print(f"<<<{final_answer}>>>")