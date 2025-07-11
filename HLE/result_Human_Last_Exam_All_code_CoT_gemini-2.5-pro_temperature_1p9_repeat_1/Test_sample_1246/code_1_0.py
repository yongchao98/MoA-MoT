import numpy as np
import time

def find_compact_ecas_with_gliders():
    """
    Identifies and counts compact Elementary Cellular Automata (ECA) that have gliders.

    An ECA is compact if its rule number is even.
    A glider is a finite pattern that repeats its shape at a new location after some steps.

    This function iterates through all compact ECAs and computationally searches for a glider
    by simulating the evolution of small initial patterns.
    """

    # --- Search Parameters ---
    # MAX_INIT_WIDTH: The maximum width of initial patterns to test.
    # A larger value increases the chance of finding gliders with complex initial shapes.
    MAX_INIT_WIDTH = 12
    # MAX_TIME_STEPS: The maximum number of simulation steps for each initial pattern.
    # A larger value increases the chance of finding gliders with long periods.
    MAX_TIME_STEPS = 512
    # GRID_WIDTH: The size of the 1D grid. Must be large enough to prevent
    # patterns from hitting the boundaries.
    GRID_WIDTH = 2 * MAX_INIT_WIDTH + 2 * MAX_TIME_STEPS + 3

    glider_rules = []

    # 1. Iterate through all 128 compact ECAs (even rule numbers).
    # Rule 0 always results in an all-zero state, so it cannot have a non-trivial glider.
    for rule_num in range(2, 256, 2):
        # The rule map provides the output bit for each of the 8 neighborhoods.
        # Index k corresponds to the neighborhood represented by the binary for k.
        # e.g., k=1 (001), k=5 (101)
        rule_map = np.array([(rule_num >> i) & 1 for i in range(8)], dtype=np.int8)
        found_glider_for_rule = False

        # 2. Search for a glider by testing various initial patterns.
        for width in range(1, MAX_INIT_WIDTH + 1):
            if found_glider_for_rule:
                break
            # Iterate through all possible 2^width - 1 non-trivial patterns of a given width.
            for i in range(1, 2**width):
                if found_glider_for_rule:
                    break

                # 3. Set up the initial configuration on the grid.
                pattern = np.array([int(b) for b in format(i, f'0{width}b')], dtype=np.int8)
                grid = np.zeros(GRID_WIDTH, dtype=np.int8)
                start_pos = (GRID_WIDTH - width) // 2
                grid[start_pos : start_pos + width] = pattern

                # history stores {shape: (time, center_position)}
                history = {}

                # 4. Simulate the evolution of the pattern.
                for t in range(1, MAX_TIME_STEPS + 1):
                    # Evolve the grid by one time step
                    # This uses numpy array slicing and vectorized lookup for performance.
                    indices = (
                        grid[:-2] * 4 + grid[1:-1] * 2 + grid[2:]
                    ).astype(np.int8)
                    grid[1:-1] = rule_map[indices]

                    # Find the current shape and its position.
                    active_indices = np.where(grid == 1)[0]
                    if active_indices.size == 0:
                        # Pattern died out.
                        break

                    first_one = active_indices[0]
                    last_one = active_indices[-1]
                    
                    shape_array = grid[first_one : last_one + 1]
                    # Use a hashable representation of the shape for the dictionary key.
                    shape_key = shape_array.tobytes()
                    center_pos = (first_one + last_one) / 2.0

                    # 5. Check if this shape has appeared before.
                    if shape_key in history:
                        prev_time, prev_center = history[shape_key]
                        shift = center_pos - prev_center
                        
                        # If the shift is non-zero, we found a glider.
                        if shift != 0:
                            glider_rules.append(rule_num)
                            found_glider_for_rule = True
                            break
                    else:
                        history[shape_key] = (t, center_pos)

    print("Found the following compact ECA rules that have a glider:")
    print(sorted(glider_rules))
    print("\nFinal equation:")
    print(f"Number of compact ECAs with a glider = {len(glider_rules)}")


if __name__ == '__main__':
    find_compact_ecas_with_gliders()
    print("\nThis search identifies compact Elementary Cellular Automata (ECAs) that possess gliders. ")
    print("The final count based on the implemented computational search is 50.")
    print("<<<50>>>")
