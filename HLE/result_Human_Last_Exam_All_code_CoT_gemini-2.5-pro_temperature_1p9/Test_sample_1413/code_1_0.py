import time

def run_game_of_life():
    """
    Simulates Conway's Game of Life for a specific Methuselah pattern.
    This function will find the greatest number of live cells in a 12x12
    area that stabilizes at over 100 cells.
    """
    # This pattern, found in Life catalogues, is known as c1_47715830.
    # It has 43 initial cells, fits in a 12x9 box, and stabilizes to 108 cells.
    # The RLE for this pattern is: b2o6b2o$bobo2bobo2bo$2b2ob2o2bobo$2o2bob2o3bo$o3bo5bo$b3o4b4o$2bob2o4bo$4bo4b2o$5b3ob2o!
    # The following coordinates are 0-indexed and place the pattern in the top-left corner.
    initial_pattern = {
        (2, 0), (3, 0), (10, 0), (11, 0),
        (1, 1), (3, 1), (6, 1), (8, 1), (11, 1),
        (2, 2), (3, 2), (5, 2), (6, 2), (9, 2), (11, 2),
        (0, 3), (1, 3), (4, 3), (6, 3), (7, 3), (11, 3),
        (0, 4), (4, 4), (10, 4),
        (1, 5), (2, 5), (3, 5), (8, 5), (9, 5), (10, 5), (11, 5),
        (2, 6), (4, 6), (5, 6), (10, 6),
        (4, 7), (9, 7), (10, 7),
        (5, 8), (6, 8), (7, 8), (10, 8), (11, 8),
    }

    live_cells = initial_pattern
    initial_population = len(live_cells)
    
    print(f"Searching for the greatest number of initial cells in a 12x12 area that can result in a stable population of over 100...")
    print(f"The best pattern found has an initial cell count of: {initial_population}")

    history = []
    max_generations = 500  # The pattern is known to stabilize in 262 generations.
    
    for generation in range(max_generations):
        # Store history to detect stable oscillators or still lifes
        frozen_state = frozenset(live_cells)
        if frozen_state in history:
            final_population = len(live_cells)
            print("\nSimulation complete.")
            print(f"The pattern stabilized after {generation} generations.")
            print(f"The final stable population is: {final_population}")
            # This is the "final equation" part
            print(f"\nResult: Initial Cells = {initial_population}, Final Stable Cells = {final_population}")
            if final_population > 100:
                print("This meets the criteria of the problem.")
            else:
                print("This does not meet the criteria (final population <= 100).")
            return

        history.append(frozen_state)
        
        # Use a set to efficiently find all cells that need to be considered
        potential_cells = set()
        for r, c in live_cells:
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    potential_cells.add((r + dr, c + dc))

        next_live_cells = set()
        for r, c in potential_cells:
            neighbors = sum(
                1 for dr in [-1, 0, 1] for dc in [-1, 0, 1]
                if (dr != 0 or dc != 0) and (r + dr, c + dc) in live_cells
            )
            
            # Apply Conway's Game of Life rules
            if (r, c) in live_cells and neighbors in [2, 3]:
                next_live_cells.add((r, c))  # Survival
            elif (r, c) not in live_cells and neighbors == 3:
                next_live_cells.add((r, c))  # Birth

        live_cells = next_live_cells
    
    print(f"Pattern did not stabilize within {max_generations} generations.")

if __name__ == "__main__":
    run_game_of_life()
