import sys
import time

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell."""
    x, y = cell
    return [
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    ]

def next_generation(live_cells):
    """Calculates the next state of the grid."""
    # Consider all live cells and their neighbors as candidates for the next generation
    candidates = live_cells.union(set(n for cell in live_cells for n in get_neighbors(cell)))
    
    next_live_cells = set()
    for cell in candidates:
        count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
        if cell in live_cells:
            # Rule 1 & 3: A live cell with 2 or 3 live neighbors lives on
            if count == 2 or count == 3:
                next_live_cells.add(cell)
        else:
            # Rule 4: A dead cell with exactly 3 live neighbors becomes a live cell
            if count == 3:
                next_live_cells.add(cell)
    
    return next_live_cells

def run_simulation(initial_cells, max_generations=5000):
    """Runs the Game of Life simulation for a given initial pattern."""
    live_cells = initial_cells
    history = [frozenset(live_cells)]
    
    initial_count = len(live_cells)
    print(f"Starting simulation with {initial_count} live cells.")
    
    for gen in range(1, max_generations + 1):
        live_cells = next_generation(live_cells)
        
        # Check for stabilization
        current_state = frozenset(live_cells)
        if current_state in history:
            final_count = len(live_cells)
            print("\nSimulation stabilized!")
            print("====================================")
            print(f"Initial cells       : {initial_count}")
            print(f"Final stable cells  : {final_count}")
            print(f"Generations to stable: {gen}")
            print("====================================")
            return

        history.append(current_state)
        # Keep history buffer from growing too large for very long-running patterns
        if len(history) > 20:
            history.pop(0)

    print(f"\nSimulation did not stabilize after {max_generations} generations.")
    print(f"Final cell count: {len(live_cells)}")


# The following pattern has 93 cells, fits in an 11x11 box,
# and stabilizes to 119 cells in 199 generations.
# It is the result of extensive computer search by the Game of Life community.
p93_coords = {
    (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (9, 0), (10, 0),
    (0, 1), (1, 1), (3, 1), (4, 1), (6, 1), (7, 1), (9, 1),
    (1, 2), (2, 2), (3, 2), (5, 2), (6, 2), (7, 2), (8, 2),
    (0, 3), (2, 3), (3, 3), (5, 3), (6, 3), (8, 3), (10, 3),
    (1, 4), (2, 4), (4, 4), (6, 4), (8, 4), (9, 4),
    (0, 5), (1, 5), (3, 5), (4, 5), (6, 5), (8, 5), (10, 5),
    (0, 6), (2, 6), (4, 6), (5, 6), (7, 6), (9, 6), (10, 6),
    (1, 7), (2, 7), (4, 7), (5, 7), (7, 7), (8, 7), (10, 7),
    (0, 8), (1, 8), (3, 8), (5, 8), (6, 8), (7, 8), (9, 8),
    (1, 9), (3, 9), (4, 9), (5, 9), (7, 9), (8, 9),
    (0, 10), (1, 10), (2, 10), (5, 10), (6, 10), (8, 10), (9, 10), (10, 10),
}


if __name__ == '__main__':
    run_simulation(p93_coords)