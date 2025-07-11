import sys

def create_pn(n):
    """Creates the initial set of live cells for the Pn pattern."""
    if n <= 0:
        return set()
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((-i, -i))
        live_cells.add((i, -i))
    return live_cells

def get_neighbors(x, y):
    """Returns the set of 8 neighbors for a given cell."""
    return {(x + dx, y + dy) for dx in [-1, 0, 1] for dy in [-1, 0, 1] if not (dx == 0 and dy == 0)}

def run_simulation_for_n(n, max_generations=250):
    """
    Runs Conway's Game of Life for a Pn pattern.

    Returns:
        A tuple (success, details) where success is a boolean and
        details is a dictionary with information about the outcome.
    """
    live_cells = create_pn(n)
    initial_dimension = 2 * n + 1
    target_dimension = 2 * initial_dimension
    history = {frozenset(live_cells)}

    for gen in range(1, max_generations + 1):
        # Determine the next generation of live cells
        candidates = live_cells.union(*(get_neighbors(x, y) for x, y in live_cells))
        next_gen_cells = set()
        for cell in candidates:
            count = len(get_neighbors(*cell) & live_cells)
            if cell in live_cells and count in (2, 3):
                next_gen_cells.add(cell)
            elif cell not in live_cells and count == 3:
                next_gen_cells.add(cell)
        
        live_cells = next_gen_cells

        # Check for termination conditions
        if not live_cells:
            return False, {"reason": "Died out", "gen": gen}

        state_snapshot = frozenset(live_cells)
        if state_snapshot in history:
            return False, {"reason": "Entered a stable cycle", "gen": gen}
        history.add(state_snapshot)

        # Check for growth condition
        min_x = min(x for x, y in live_cells)
        max_x = max(x for x, y in live_cells)
        min_y = min(y for x, y in live_cells)
        max_y = max(y for x, y in live_cells)

        current_width = max_x - min_x + 1
        current_height = max_y - min_y + 1

        if current_width >= target_dimension or current_height >= target_dimension:
            return True, {
                "gen": gen,
                "n": n,
                "initial_dim": initial_dimension,
                "target_dim": target_dimension,
                "final_width": current_width,
                "final_height": current_height
            }

    return False, {"reason": f"Exceeded {max_generations} generations", "gen": max_generations}


def find_smallest_growing_pn():
    """Finds the smallest n > 0 for which Pn grows to twice its size."""
    n = 1
    while True:
        sys.stdout.write(f"Testing n={n}...")
        sys.stdout.flush()
        
        success, details = run_simulation_for_n(n)
        
        if success:
            print(f"\nSuccess! Found the smallest value n > 0 is {details['n']}.")
            
            initial_dim = details['initial_dim']
            target_dim = details['target_dim']
            final_dim = max(details['final_width'], details['final_height'])
            
            print(f"The P{n} pattern has an initial dimension of {initial_dim}.")
            print(f"The target dimension was >= 2 * {initial_dim}, which is {target_dim}.")
            print(f"At generation {details['gen']}, the pattern's dimension grew to {final_dim}.")
            print("\nFinal equation check:")
            # Output each number in the final equation
            print(f"{final_dim} >= 2 * (2 * {n} + 1)")
            print(f"{final_dim} >= 2 * ({initial_dim})")
            print(f"{final_dim} >= {target_dim}")
            print(f"This is True, so n={n} is the answer.")
            return n
        else:
            print(f" Failed ({details['reason']})")
            n += 1

if __name__ == '__main__':
    find_smallest_growing_pn()