import sys

def get_neighbors(cell):
    """Returns a set of 8 neighbor coordinates for a given cell."""
    x, y = cell
    return {
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    }

def step(live_cells):
    """Computes the next generation of live cells based on Conway's rules."""
    # The set of candidates for the next generation includes all live cells
    # and all their immediate neighbors.
    candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    
    next_gen_live_cells = set()
    
    for cell in candidates:
        live_neighbors_count = len(get_neighbors(cell).intersection(live_cells))
        
        # Rule for a currently live cell
        if cell in live_cells:
            if live_neighbors_count == 2 or live_neighbors_count == 3:
                next_gen_live_cells.add(cell)
        # Rule for a currently dead cell (reproduction)
        else:
            if live_neighbors_count == 3:
                next_gen_live_cells.add(cell)
                
    return next_gen_live_cells

def get_dimensions(live_cells):
    """Calculates the width and height of the bounding box around live cells."""
    if not live_cells:
        return 0, 0
    min_x = min(x for x, y in live_cells)
    max_x = max(x for x, y in live_cells)
    min_y = min(y for x, y in live_cells)
    max_y = max(y for x, y in live_cells)
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    return width, height

def create_pn_pattern(n):
    """Creates the initial set of live cells for the Pn pattern."""
    if n <= 0:
        return set()
    live_cells = {(0, 0)}
    for i in range(1, n + 1):
        live_cells.add((i, i))
        live_cells.add((-i, i))
        live_cells.add((i, -i))
        live_cells.add((-i, -i))
    return live_cells

def solve():
    """
    Finds the smallest integer n > 0 for which the Pn pattern grows
    to at least twice its original size in any dimension.
    """
    n = 1
    # A safety limit to prevent infinite loops for chaotic patterns
    max_generations = 500
    
    while True:
        live_cells = create_pn_pattern(n)
        initial_size = 2 * n + 1
        target_size = 2 * initial_size

        history = [live_cells]
        
        for gen in range(1, max_generations + 1):
            live_cells = step(live_cells)
            
            # Condition 1: Pattern dies out
            if not live_cells:
                break
            # Condition 2: Pattern stabilizes or enters a cycle
            if live_cells in history:
                break
            
            history.append(live_cells)
            
            width, height = get_dimensions(live_cells)
            
            # Condition 3: Growth condition is met
            if width >= target_size or height >= target_size:
                final_dimension = max(width, height)
                
                print(f"The smallest value for n is {n}.")
                print("\nThis was found by checking the condition: final_dimension >= 2 * (2 * n + 1)")
                print(f"For n={n}, the initial size is {initial_size} and the target size is {target_size}.")
                print(f"After {gen} generations, the pattern grew to a size of {width}x{height}.")
                print(f"Substituting the values into the check: {final_dimension} >= 2 * (2 * {n} + 1)")

                print("\nThe numbers in the final verification equation are:")
                print(final_dimension)
                print(2)
                print(2)
                print(n)
                print(1)

                return n
        
        n += 1
        # Safety break for the search loop
        if n > 20:
            print("Could not find a solution for n up to 20.")
            return None

# Execute the solver function to find and print the result.
result = solve()
if result is not None:
    print(f"\n<<< {result} >>>")
