import collections

def find_smallest_growing_pn():
    """
    Finds the smallest integer n > 0 for which the Pn pattern in Conway's Game of Life
    grows to at least twice its original size in any dimension.
    """

    def generate_pn(n):
        """Generates the initial live cells for a Pn pattern."""
        live_cells = set()
        live_cells.add((0, 0))
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, -i))
            live_cells.add((i, -i))
            live_cells.add((-i, i))
        return live_cells

    def evolve(live_cells):
        """Computes the next generation of live cells based on the rules of the game."""
        if not live_cells:
            return set()

        # Use a Counter to store the number of live neighbors for each cell
        neighbor_counts = collections.Counter()
        for r, c in live_cells:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    neighbor_counts[(r + dr, c + dc)] += 1

        # Apply Game of Life rules to determine the next generation of live cells
        next_live_cells = set()
        for cell, count in neighbor_counts.items():
            # A cell becomes live if it has exactly 3 neighbors (reproduction)
            if count == 3:
                next_live_cells.add(cell)
            # A cell stays live if it has 2 neighbors and was already live (survival)
            elif count == 2 and cell in live_cells:
                next_live_cells.add(cell)
        return next_live_cells

    def get_bounding_box_size(live_cells):
        """Calculates the width and height of the bounding box of live cells."""
        if not live_cells:
            return 0, 0
        min_r = min(cell[0] for cell in live_cells)
        max_r = max(cell[0] for cell in live_cells)
        min_c = min(cell[1] for cell in live_cells)
        max_c = max(cell[1] for cell in live_cells)
        # Width is based on columns (x-axis), Height is based on rows (y-axis)
        width = max_c - min_c + 1
        height = max_r - min_r + 1
        return width, height

    # Iterate through n starting from 1 to find the smallest solution
    n = 0
    while True:
        n += 1
        live_cells = generate_pn(n)
        
        # Calculate initial and target sizes
        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        # Keep a history to detect stable or oscillating patterns
        history = {frozenset(live_cells)}
        
        # Limit generations to prevent truly infinite computations,
        # though history check should handle non-growing patterns.
        max_generations = 1000

        for g in range(1, max_generations + 1):
            live_cells = evolve(live_cells)
            
            # Condition 1: Pattern died out
            if not live_cells:
                break
            
            # Condition 2: Pattern is stable or oscillating
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)
            
            # Condition 3: Check if the pattern has grown to the target size
            width, height = get_bounding_box_size(live_cells)
            if width >= target_size or height >= target_size:
                print(f"Found solution for n = {n}")
                print(f"The initial pattern P{n} has a size of {initial_size}x{initial_size}.")
                print(f"The equation for the growth condition is: width >= 2 * {initial_size} or height >= 2 * {initial_size}")
                print(f"The target size is {target_size}.")
                print(f"After {g} generations, the pattern's bounding box grew to {width}x{height}.")
                print(f"Final equation check: {width} >= {target_size} or {height} >= {target_size} is True.")
                return n
                
# Run the search and print the final result
solution_n = find_smallest_growing_pn()
print(f"\nThe smallest value for n is {solution_n}.")
<<<5>>>