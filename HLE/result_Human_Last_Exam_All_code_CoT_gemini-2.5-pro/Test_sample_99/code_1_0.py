import collections

def get_bounding_box(live_cells):
    """Calculates the bounding box width and height for a set of live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def run_generation(live_cells):
    """Computes the next generation of live cells based on Conway's rules."""
    neighbor_counts = collections.Counter()
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                neighbor_counts[(x + dx, y + dy)] += 1

    next_live_cells = set()
    for cell, count in neighbor_counts.items():
        if count == 3 and cell not in live_cells:
            next_live_cells.add(cell)
        elif count in [2, 3] and cell in live_cells:
            next_live_cells.add(cell)
    return next_live_cells

def find_smallest_growing_pn():
    """
    Finds the smallest integer n > 0 for pattern Pn that grows to at least
    twice its original size in any dimension.
    """
    max_n_to_check = 10
    max_generations = 200

    for n in range(1, max_n_to_check + 1):
        # 1. Define the initial Pn pattern
        live_cells = set([(0, 0)])
        for i in range(1, n + 1):
            live_cells.add((i, i))
            live_cells.add((-i, i))
            live_cells.add((i, -i))
            live_cells.add((-i, -i))
        
        # 2. Calculate initial and target sizes
        initial_size = 2 * n + 1
        target_size = 2 * initial_size
        
        current_cells = live_cells
        history = {frozenset(current_cells)}
        
        # 3. Run the simulation
        for gen in range(1, max_generations + 1):
            current_cells = run_generation(current_cells)
            
            if not current_cells:
                break
                
            frozen_state = frozenset(current_cells)
            if frozen_state in history:
                break
            history.add(frozen_state)

            # 4. Check the growth condition
            width, height = get_bounding_box(current_cells)
            if width >= target_size or height >= target_size:
                print(f"The smallest value for n is {n}.")
                print(f"The initial pattern P{n} had an original size of {initial_size}x{initial_size}.")
                print(f"The target size was >= {target_size} in any dimension.")
                print(f"After {gen} generations, the pattern grew to a size of {width}x{height}, meeting the condition.")
                return n
                
    print(f"No solution found for n up to {max_n_to_check} within {max_generations} generations.")
    return None

if __name__ == '__main__':
    find_smallest_growing_pn()