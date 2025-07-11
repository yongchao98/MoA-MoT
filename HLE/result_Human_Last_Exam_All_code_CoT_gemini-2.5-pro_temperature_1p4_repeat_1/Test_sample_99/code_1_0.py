import collections

def solve_game_of_life_pn():
    """
    Finds the smallest integer n > 0 for Conway's Game of Life pattern Pn 
    such that the pattern grows to at least twice its original size in any dimension.
    """

    def get_p_n_pattern(n):
        """Generates the initial set of live cells for the Pn pattern."""
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))
        return live_cells

    def get_bounding_box_size(cells):
        """Calculates the width and height of the pattern's bounding box."""
        if not cells:
            return 0, 0
        min_x = min(c[0] for c in cells)
        max_x = max(c[0] for c in cells)
        min_y = min(c[1] for c in cells)
        max_y = max(c[1] for c in cells)
        width = max_x - min_x + 1
        height = max_y - min_y + 1
        return width, height

    def run_one_generation(live_cells):
        """Simulates one generation (step) of Conway's Game of Life."""
        if not live_cells:
            return set()
        
        # Count live neighbors for all relevant cells
        neighbor_counts = collections.defaultdict(int)
        for x, y in live_cells:
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    neighbor_counts[(x + dx, y + dy)] += 1

        next_gen_cells = set()
        
        # Apply Game of Life rules
        # A cell is alive in the next generation if:
        # - It is currently alive and has 2 or 3 neighbors.
        # - It is currently dead and has exactly 3 neighbors.
        # This simplifies to: any cell with 3 neighbors is alive, or any
        # currently live cell with 2 neighbors stays alive.
        
        cells_to_check = set(neighbor_counts.keys())
        for cell in cells_to_check:
            count = neighbor_counts[cell]
            if count == 3:
                next_gen_cells.add(cell)
            elif count == 2 and cell in live_cells:
                next_gen_cells.add(cell)
                
        return next_gen_cells

    n = 1
    max_generations_to_check = 500  # A safety limit to prevent very long runs

    while True:
        initial_cells = get_p_n_pattern(n)
        initial_size, _ = get_bounding_box_size(initial_cells) # Width and height are the same
        target_size = 2 * initial_size
        
        current_cells = initial_cells
        history = {frozenset(current_cells)}

        for gen in range(1, max_generations_to_check + 1):
            current_cells = run_one_generation(current_cells)
            
            # Condition 1: Pattern died out
            if not current_cells:
                break

            # Condition 2: Pattern is stable or oscillating
            frozen_state = frozenset(current_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            width, height = get_bounding_box_size(current_cells)
            
            # Condition 3: Growth target reached
            if width >= target_size or height >= target_size:
                print(f"Solution found for n = {n}")
                print(f"The initial size of P{n} is {initial_size}x{initial_size}.")
                print(f"The target size (2 * initial size) is {target_size}.")
                print(f"At generation {gen}, the pattern reached a size of {width}x{height}.")
                print(f"\nThe smallest value for n is {n}.")
                return n

        n += 1

if __name__ == '__main__':
    solve_game_of_life_pn()