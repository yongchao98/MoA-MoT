import collections

def create_pn(n):
    """Creates the set of live cells for a Pn pattern."""
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def get_bounding_box_and_size(live_cells):
    """Calculates the bounding box and size of a set of live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(c[0] for c in live_cells)
    max_x = max(c[0] for c in live_cells)
    min_y = min(c[1] for c in live_cells)
    max_y = max(c[1] for c in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def step(live_cells):
    """Computes the next generation of live cells."""
    if not live_cells:
        return set()

    # This counter will store potential new cells and the number of live neighbors they have.
    neighbor_counts = collections.Counter()
    for x, y in live_cells:
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                if dx == 0 and dy == 0:
                    continue
                neighbor_counts[(x + dx, y + dy)] += 1
    
    next_generation = set()
    
    # Iterate through all cells that have at least one live neighbor
    for cell, count in neighbor_counts.items():
        # Rule: A dead cell with exactly three live neighbours becomes a live cell.
        if count == 3 and cell not in live_cells:
            next_generation.add(cell)
            
    # Iterate through the original live cells
    for cell in live_cells:
        # Rule: A live cell with two or three live neighbours lives on.
        if neighbor_counts.get(cell, 0) in [2, 3]:
            next_generation.add(cell)
            
    return next_generation

def find_smallest_growing_pn():
    """
    Finds the smallest integer n > 0 for a Pn pattern that grows to
    at least twice its original size.
    """
    n = 1
    max_generations = 500  # A safe limit for generations per n

    while True:
        live_cells = create_pn(n)
        
        # Initial and target sizes
        original_dimension = 2 * n + 1
        target_dimension = 2 * original_dimension
        
        history = {frozenset(live_cells)}

        for gen in range(1, max_generations + 1):
            live_cells = step(live_cells)

            # Condition 1: Pattern died out
            if not live_cells:
                break

            # Condition 2: Pattern grows to target size
            current_width, current_height = get_bounding_box_and_size(live_cells)
            
            growth_found = False
            satisfied_dimension = 0
            if current_width >= target_dimension:
                growth_found = True
                satisfied_dimension = current_width
            elif current_height >= target_dimension:
                growth_found = True
                satisfied_dimension = current_height

            if growth_found:
                print(f"Solution found for n = {n}")
                print(f"The pattern grew to at least twice its original size at generation {gen}.")
                print(f"Original dimension: {original_dimension}")
                print(f"Target dimension: {target_dimension}")
                print("The condition was met as follows:")
                # Print the final equation with each number
                print(f"{satisfied_dimension} >= 2 * (2 * {n} + 1)")
                return n

            # Condition 3: Pattern stabilized or is an oscillator
            frozen_state = frozenset(live_cells)
            if frozen_state in history:
                break 
            history.add(frozen_state)
        
        n += 1

# Execute the function to find the answer
solution_n = find_smallest_growing_pn()
print(f"\nThe smallest value for Pn > 0 is {solution_n}.")
print(f"<<<{solution_n}>>>")