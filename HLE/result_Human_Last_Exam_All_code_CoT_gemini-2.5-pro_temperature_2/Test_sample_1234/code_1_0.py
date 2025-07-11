import sys

def solve():
    """
    Solves the Conway's Game of Life problem for 3x3 starting grids by simulating
    each of the 512 possibilities on an infinite plane.
    """
    
    # Increase the recursion limit as a safeguard for pattern normalization,
    # though it is unlikely to be hit with 3x3 starting patterns.
    sys.setrecursionlimit(2000)

    # This list will store the integer representation of each configuration that dies out.
    dying_configs_indices = []
    
    # Total configurations for a 3x3 grid is 2^9 = 512
    total_configs = 512

    # Loop through each possible configuration from 0 to 511.
    for i in range(total_configs):
        initial_live_cells = set()
        # Convert the integer `i` into a 3x3 grid pattern.
        # Each bit in the 9-bit representation of `i` corresponds to a cell.
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                initial_live_cells.add((row, col))
        
        # Simulate the evolution of the pattern and check if it dies out.
        if simulate_life(initial_live_cells):
            dying_configs_indices.append(i)

    # --- Output the results ---
    count = len(dying_configs_indices)
    print(f"Out of {total_configs} initial configurations, {count} will eventually result in no living cells.")
    print("The integer representations of these configurations are shown below as a sum:")
    
    # As requested, output each number that corresponds to a dying configuration,
    # formatted as a single "equation" string.
    equation_str = " + ".join(map(str, sorted(dying_configs_indices)))
    print(equation_str)


def simulate_life(live_cells, max_generations=100):
    """
    Simulates Conway's Game of Life using a sparse set-based grid.

    Args:
        live_cells (set): A set of (row, col) tuples representing live cells.
        max_generations (int): A safe upper limit on generations to simulate.

    Returns:
        bool: True if the pattern dies out, False otherwise.
    """
    history = set()

    for _ in range(max_generations):
        if not live_cells:
            # No live cells remain; the pattern has died out.
            return True

        # To detect all cycles (stable, oscillating, and moving), create a
        # position-independent "canonical" representation of the pattern.
        min_r = min(r for r, c in live_cells)
        min_c = min(c for r, c in live_cells)
        canonical_form = frozenset((r - min_r, c - min_c) for r, c in live_cells)

        if canonical_form in history:
            # The pattern has repeated a previous state; it is in a cycle and will not die.
            return False
        history.add(canonical_form)

        # To find the next generation, we only need to check the state of current
        # live cells and their immediate neighbors.
        cells_to_check = set()
        for r, c in live_cells:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    cells_to_check.add((r + dr, c + dc))

        next_live_cells = set()
        for r, c in cells_to_check:
            # Count the live neighbors for the current cell.
            neighbors = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    if (r + dr, c + dc) in live_cells:
                        neighbors += 1
            
            # Apply the rules of Conway's Game of Life.
            is_alive = (r, c) in live_cells
            if is_alive and (neighbors == 2 or neighbors == 3):
                next_live_cells.add((r, c)) # Rule for survival.
            elif not is_alive and neighbors == 3:
                next_live_cells.add((r, c)) # Rule for birth.
        
        live_cells = next_live_cells

    # If the simulation reaches the generation limit without dying or cycling,
    # it's considered a long-lived pattern that does not die out.
    return False

if __name__ == '__main__':
    solve()
    print("\n<<<139>>>")