import collections

def find_smallest_growing_pn():
    """
    Solves for the smallest integer n > 0 for which the Conway's Game of Life
    pattern Pn grows to at least twice its original size in any dimension.
    """
    n = 1
    # A safety limit on generations to prevent any unforeseen infinite loops.
    # Cycle detection should handle most cases.
    max_generations = 2000

    while True:
        # Step 1: Define the Pn pattern and its initial/target sizes.
        # A set of (x,y) tuples represents the live cells.
        live_cells = set()
        if n > 0:
            live_cells.add((0, 0))
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))

        # The initial bounding box for Pn is a square.
        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        history = {frozenset(live_cells)}
        
        # Step 2: Simulate generations for the current value of n.
        for generation in range(1, max_generations + 1):
            
            # Count neighbors for all relevant cells.
            neighbor_counts = collections.defaultdict(int)
            for x, y in live_cells:
                for i in range(-1, 2):
                    for j in range(-1, 2):
                        if i == 0 and j == 0:
                            continue
                        neighbor_counts[(x + i, y + j)] += 1

            # Apply Game of Life rules to get the next state.
            next_live_cells = set()
            for cell, count in neighbor_counts.items():
                # A cell is alive in the next generation if:
                # 1. It was dead and has exactly 3 neighbors (reproduction).
                # 2. It was alive and has 2 or 3 neighbors (survival).
                if count == 3 or (count == 2 and cell in live_cells):
                    next_live_cells.add(cell)
            
            live_cells = next_live_cells

            # Step 3: Check for termination or success conditions.
            if not live_cells:
                # The pattern has died out.
                break

            # Check if the pattern has become stable or periodic.
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            # Get the current size of the pattern's bounding box.
            min_x = min(x for x, y in live_cells)
            max_x = max(x for x, y in live_cells)
            min_y = min(y for x, y in live_cells)
            max_y = max(y for x, y in live_cells)
            current_width = max_x - min_x + 1
            current_height = max_y - min_y + 1

            # Check if the growth condition is met.
            if current_width >= target_size or current_height >= target_size:
                print(f"Solution found for Pn with n = {n}.")
                print(f"The initial size of the pattern was {initial_size}x{initial_size}.")
                print(f"The target size was a width or height of at least {target_size}.")
                
                print("\nHere are the numbers in the equations used for the final check:")
                print(f"Initial size = (2 * n) + 1 = (2 * {n}) + 1 = {initial_size}")
                print(f"Target size = 2 * Initial size = 2 * {initial_size} = {target_size}")
                print(f"At generation {generation}, the bounding box grew to {current_width}x{current_height}.")
                print(f"This met the condition because {max(current_width, current_height)} >= {target_size}.")
                print(f"\nTherefore, the smallest value for n > 0 is {n}.")
                return # Exit after finding the solution.
        
        # Step 4: If the simulation for the current n ended without success, try the next n.
        n += 1

if __name__ == '__main__':
    find_smallest_growing_pn()