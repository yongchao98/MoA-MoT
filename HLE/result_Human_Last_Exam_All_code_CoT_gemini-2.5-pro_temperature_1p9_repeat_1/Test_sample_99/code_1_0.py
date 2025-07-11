import sys

def create_pn(n):
    """Creates the initial set of live cells for a Pn pattern."""
    if n <= 0:
        return set()
    cells = {(0, 0)}
    for i in range(1, n + 1):
        cells.add((i, i))
        cells.add((-i, i))
        cells.add((i, -i))
        cells.add((-i, -i))
    return cells

def get_neighbors(cell):
    """Returns the 8 neighbors of a given cell."""
    x, y = cell
    return {
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    }

def evolve(live_cells):
    """Calculates the next generation of live cells."""
    # Candidates for the next generation are the current live cells
    # and their immediate neighbors.
    candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    next_live_cells = set()
    for cell in candidates:
        count = len(get_neighbors(cell).intersection(live_cells))
        is_live = cell in live_cells
        
        # Rule 1 & 2: A live cell with 2 or 3 neighbors survives.
        if is_live and count in [2, 3]:
            next_live_cells.add(cell)
        # Rule 4: A dead cell with exactly 3 neighbors becomes a live cell.
        elif not is_live and count == 3:
            next_live_cells.add(cell)
            
    return next_live_cells

def get_bounding_box(live_cells):
    """Calculates the width and height of the pattern."""
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def solve_game_of_life_problem():
    """
    Iterates through n > 0 to find the smallest Pn pattern that grows
    to at least twice its original size.
    """
    n = 0
    # A limit to prevent true infinite loops, though history tracking should handle it.
    max_n_to_check = 20 

    while n < max_n_to_check:
        n += 1
        live_cells = create_pn(n)
        
        initial_dim = 2 * n + 1
        target_dim = 2 * initial_dim

        # Keep a history of states to detect stable or oscillating patterns
        history = [frozenset(live_cells)]
        
        # Max generations to check for a single n
        for generation in range(1, 500):
            live_cells = evolve(live_cells)
            
            # Check for termination conditions
            if not live_cells:
                break  # Pattern died out
                
            current_state = frozenset(live_cells)
            if current_state in history:
                break # Pattern is stable or oscillating
            history.append(current_state)

            # Check for growth condition
            width, height = get_bounding_box(live_cells)
            if width >= target_dim or height >= target_dim:
                print(f"Success! Found the smallest value for Pn > 0.")
                print(f"The answer is for the pattern Pn where n = {n}.")
                print(f"The initial dimensions of P{n} were {initial_dim}x{initial_dim}.")
                print(f"The required dimension was >= {target_dim}.")
                print(f"At generation {generation}, the pattern grew to {width}x{height}.")
                print("\nThe smallest value for n is:")
                print(n)
                return

    print("No solution found within the checked range.")

if __name__ == '__main__':
    solve_game_of_life_problem()