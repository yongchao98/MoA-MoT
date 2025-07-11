import sys

def find_smallest_growing_pn():
    """
    Solves for the smallest integer n > 0 where the Game of Life pattern Pn
    grows to at least twice its initial size along any dimension.

    A Pn pattern is an X-shape with a central cell and n live cells
    radiating diagonally from each of its four corners.
    """
    n = 1
    # We set a practical limit for the search, though the answer is found quickly.
    max_n_to_check = 100

    while n <= max_n_to_check:
        # Step 1: Generate the initial live cells for the Pn pattern.
        live_cells = {(0, 0)}
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))

        # Step 2: Calculate the initial and target dimensions.
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # Keep track of recent states to detect stable or oscillating patterns.
        previous_states = []
        max_generations = 500  # Max simulation steps for each n.

        for generation in range(max_generations):
            # Check for extinction (pattern dies out).
            if not live_cells:
                break  # This n fails, move to the next.

            # Use a frozenset for efficient checking of previous states.
            current_state_fs = frozenset(live_cells)
            if current_state_fs in previous_states:
                break  # Pattern is stable or oscillating, so this n fails.
            
            # Store the last 10 states to detect common oscillators.
            previous_states.append(current_state_fs)
            if len(previous_states) > 10:
                previous_states.pop(0)

            # Step 3: Check if the growth condition has been met.
            min_x = min(cell[0] for cell in live_cells)
            max_x = max(cell[0] for cell in live_cells)
            min_y = min(cell[1] for cell in live_cells)
            max_y = max(cell[1] for cell in live_cells)

            current_width = max_x - min_x + 1
            current_height = max_y - min_y + 1

            if current_width >= target_dim or current_height >= target_dim:
                print(f"Success! The smallest value for n is {n}.")
                print(f"The initial P{n} pattern had dimensions {initial_dim}x{initial_dim}.")
                print(f"The growth target was a dimension of at least {target_dim}.")
                print(f"At generation {generation + 1}, the pattern reached a size of {current_width}x{current_height}.")
                
                # Showing the final check as an "equation"
                if current_width >= target_dim:
                    print(f"Final check: Width ({current_width}) >= Target ({target_dim})")
                if current_height >= target_dim:
                    print(f"Final check: Height ({current_height}) >= Target ({target_dim})")
                return n

            # Step 4: Calculate the next generation of live cells.
            # We only need to check cells that are live or are neighbors of live cells.
            candidates_to_check = set()
            for r, c in live_cells:
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        candidates_to_check.add((r + dr, c + dc))
            
            next_live_cells = set()
            for r, c in candidates_to_check:
                live_neighbors = 0
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        if dr == 0 and dc == 0:
                            continue
                        if (r + dr, c + dc) in live_cells:
                            live_neighbors += 1
                
                # Apply the rules of Life
                is_currently_live = (r, c) in live_cells
                if is_currently_live and live_neighbors in [2, 3]:
                    next_live_cells.add((r, c))  # Survival
                elif not is_currently_live and live_neighbors == 3:
                    next_live_cells.add((r, c))  # Birth
            
            live_cells = next_live_cells
        
        # If the simulation loop completes without success, try the next n.
        n += 1

    print(f"No solution was found for n up to {max_n_to_check}.")
    return None

if __name__ == '__main__':
    find_smallest_growing_pn()