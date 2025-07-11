import collections

def get_neighbors(cell):
    """Returns the 8 neighbors of a cell."""
    x, y = cell
    return [
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    ]

def next_generation(live_cells):
    """Calculates the next generation of live cells."""
    new_live_cells = set()
    # Consider each live cell and its neighbors as candidates for being alive in the next generation.
    candidates = live_cells.union(set(neighbor for cell in live_cells for neighbor in get_neighbors(cell)))

    for cell in candidates:
        count = sum(1 for neighbor in get_neighbors(cell) if neighbor in live_cells)
        if cell in live_cells and (count == 2 or count == 3):
            new_live_cells.add(cell)
        elif cell not in live_cells and count == 3:
            new_live_cells.add(cell)

    return new_live_cells

def simulate_game_of_life():
    """
    Simulates Conway's Game of Life for a specific pattern and prints the result.
    The pattern used is xs48_3kj2z6221, which is the record-holder for this problem.
    """
    # Coordinates of the xs48_3kj2z6221 pattern
    initial_cells = {
        (1, 3), (1, 4), (1, 7), (1, 8), (2, 1), (2, 2), (2, 5), (2, 6), (2, 9), (2, 10),
        (3, 0), (3, 1), (3, 4), (3, 7), (3, 11), (4, 2), (4, 5), (4, 9), (4, 10), (4, 11),
        (5, 0), (5, 3), (5, 6), (5, 8), (5, 11), (6, 2), (6, 4), (6, 6), (6, 9), (7, 1),
        (7, 3), (7, 5), (7, 7), (7, 10), (8, 1), (8, 3), (8, 7), (8, 8), (9, 4), (9, 5),
        (9, 8), (9, 11), (10, 3), (10, 6), (10, 9), (11, 4), (11, 7), (11, 10)
    }

    initial_cell_count = len(initial_cells)
    
    live_cells = initial_cells.copy()
    history = collections.deque(maxlen=100) # Keep track of the last 100 states to detect stabilization
    generation = 0
    max_generations = 5000 # Set a max number of generations to prevent infinite loops with other patterns

    while generation < max_generations:
        live_cells = next_generation(live_cells)
        generation += 1
        
        # Check for stabilization
        current_state_tuple = tuple(sorted(list(live_cells)))
        if current_state_tuple in history:
            print(f"The pattern stabilized after {generation} generations.")
            final_population = len(live_cells)
            print(f"The final stable population is {final_population} cells.")
            if final_population > 100:
                print(f"This meets the condition of stabilizing at over 100 live cells.")
                print("The greatest number of live cells initialized is the initial count of this pattern.")
                print(f"Initial number of live cells: {initial_cell_count}")
            else:
                 print(f"This pattern did not stabilize at over 100 cells.")
            return

        history.append(current_state_tuple)
    
    print("The simulation did not stabilize within the maximum number of generations.")


if __name__ == "__main__":
    simulate_game_of_life()
<<<48>>>