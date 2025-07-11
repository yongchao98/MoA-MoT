import sys

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell."""
    x, y = cell
    return {
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    }

def step(live_cells):
    """Computes the next generation of live cells."""
    if not live_cells:
        return set()

    # We only need to check the live cells and their immediate neighbors
    candidates = live_cells.union(*(get_neighbors(cell) for cell in live_cells))
    
    next_generation = set()
    for cell in candidates:
        count = len(get_neighbors(cell).intersection(live_cells))
        
        # Apply Conway's Game of Life rules
        if cell in live_cells and count in {2, 3}:
            next_generation.add(cell)
        elif cell not in live_cells and count == 3:
            next_generation.add(cell)
            
    return next_generation

def run_simulation(initial_pattern, max_generations=500):
    """
    Runs the Game of Life simulation for a given pattern.
    
    The function prints the initial and final state and returns the counts.
    """
    # The "Century" pattern has 11 cells and stabilizes to 101 cells in 100 generations.
    # It fits within a 7x5 box.
    # RLE: 3b2o$3bo2bo$2obo$bo2bobo$b2o!
    pattern_coords = {
        (3, 0), (4, 0),
        (3, 1), (5, 1),
        (2, 2), (4, 2),
        (1, 3), (4, 3), (6, 3),
        (1, 4), (2, 4),
    }

    live_cells = pattern_coords
    initial_count = len(live_cells)
    
    print(f"Simulating the '{initial_pattern}' pattern.")
    print(f"Initial live cells: {initial_count}")

    for i in range(max_generations):
        next_cells = step(live_cells)
        
        # Check for stabilization (pattern becomes a still life)
        if next_cells == live_cells:
            final_count = len(live_cells)
            print(f"Stabilized after {i} generations.")
            print(f"Final live cells: {final_count}")
            print("---")
            print(f"Resulting equation: {initial_count} -> {final_count}")
            return
            
        live_cells = next_cells

    print(f"Did not stabilize within {max_generations} generations.")
    final_count = len(live_cells)
    print(f"Live cells at end: {final_count}")

# Run the simulation for the "Century" pattern
run_simulation("Century")
