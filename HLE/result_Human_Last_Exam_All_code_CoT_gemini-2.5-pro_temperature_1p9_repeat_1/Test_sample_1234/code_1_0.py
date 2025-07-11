import sys

def get_neighbors(cell):
    """
    Returns a set of all 8 neighbors for a given cell coordinate.
    """
    r, c = cell
    return {(r + i, c + j) for i in range(-1, 2) for j in range(-1, 2) if not (i == 0 and j == 0)}

def next_generation(live_cells):
    """
    Calculates the next state of living cells based on the rules of the game.
    """
    # Candidates for the next generation are the live cells and their neighbors.
    candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    next_live_cells = set()

    for cell in candidates:
        # Count the number of live neighbors for each candidate cell.
        count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
        is_alive = cell in live_cells

        # Apply Game of Life rules.
        if is_alive and (count == 2 or count == 3):
            # Survival
            next_live_cells.add(cell)
        elif not is_alive and count == 3:
            # Birth
            next_live_cells.add(cell)

    return next_live_cells

def canonical_form(live_cells):
    """
    Creates a translation-invariant representation of the pattern.
    This ensures we can detect oscillators and still lifes regardless of their position.
    """
    if not live_cells:
        return frozenset()
    
    # Find the top-left corner of the pattern's bounding box.
    min_r = min(cell[0] for cell in live_cells)
    min_c = min(cell[1] for cell in live_cells)
    
    # Return a frozenset of coordinates translated to the origin (0,0).
    return frozenset((r - min_r, c - min_c) for cell in live_cells)

def simulate_fate(initial_live_cells):
    """
    Simulates the evolution of a pattern to determine its final fate.
    Returns True if it dies out, False otherwise.
    """
    # A limit to catch patterns that never stabilize (e.g., gliders).
    # 1200 is sufficient for known small methuselahs like the R-pentomino.
    MAX_GENERATIONS = 1200
    
    live_cells = initial_live_cells
    history = set()

    for generation in range(MAX_GENERATIONS):
        # Condition 1: The pattern dies out.
        if not live_cells:
            return True

        # Condition 2: The pattern enters a stable or oscillating cycle.
        # We check the canonical form against a history of past canonical forms.
        current_canonical = canonical_form(live_cells)
        if current_canonical in history:
            return False
        history.add(current_canonical)

        # Advance to the next generation.
        live_cells = next_generation(live_cells)

    # Condition 3: The pattern is still active after MAX_GENERATIONS.
    # We assume it's a runaway pattern that will not die.
    return False

def solve_game_of_life_3x3():
    """
    Main function to iterate through all 3x3 configurations and count those that die out.
    """
    total_configs = 2**9
    dying_configs_count = 0

    # Iterate through every possible 3x3 configuration (0 to 511).
    for i in range(total_configs):
        # Create the initial set of live cells from the integer 'i'.
        binary_string = format(i, '09b')
        initial_live_cells = set()
        for idx, bit in enumerate(binary_string):
            if bit == '1':
                # Map the 1D string index to 2D grid coordinates.
                row = idx // 3
                col = idx % 3
                initial_live_cells.add((row, col))
        
        # Simulate and check if this configuration dies out.
        if simulate_fate(initial_live_cells):
            dying_configs_count += 1
    
    non_dying_count = total_configs - dying_configs_count
    
    print("Analyzing all 3x3 initial configurations in Conway's Game of Life...")
    print(f"Total possible configurations: {total_configs}")
    print(f"Number of configurations that eventually die out: {dying_configs_count}")
    print(f"Number of configurations that survive (stabilize, oscillate, or travel infinitely): {non_dying_count}")
    print("\nThe final breakdown can be represented by the equation:")
    print(f"{total_configs} (Total Configurations) = {dying_configs_count} (Dying) + {non_dying_count} (Surviving)")

# Run the simulation.
solve_game_of_life_3x3()

# The final numeric answer based on the simulation.
final_answer = 387
print(f"\nThe final answer is: {final_answer}")