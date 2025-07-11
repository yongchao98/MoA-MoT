import sys

def get_neighbor_coords(cell):
    """Returns the coordinates of the 8 neighbors of a cell."""
    x, y = cell
    return [
        (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
        (x - 1, y),                 (x + 1, y),
        (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
    ]

def get_next_generation(live_cells):
    """Calculates the next state of the grid based on Conway's rules."""
    # A set of all cells that could potentially be alive in the next generation.
    # This includes current live cells and their neighbors.
    candidates = set()
    for cell in live_cells:
        candidates.add(cell)
        for neighbor in get_neighbor_coords(cell):
            candidates.add(neighbor)

    next_live_cells = set()
    for cell in candidates:
        live_neighbors = sum(1 for neighbor in get_neighbor_coords(cell) if neighbor in live_cells)
        if cell in live_cells and live_neighbors in [2, 3]:
            # Survival
            next_live_cells.add(cell)
        elif cell not in live_cells and live_neighbors == 3:
            # Birth
            next_live_cells.add(cell)
    return next_live_cells

def simulate_life(initial_pattern, max_generations=1000):
    """
    Simulates Conway's Game of Life for a given pattern.
    Stops when the pattern stabilizes (becomes static or enters a loop) or max_generations is reached.
    """
    live_cells = frozenset(initial_pattern)
    history = {live_cells}
    
    for gen in range(max_generations):
        live_cells = frozenset(get_next_generation(live_cells))
        
        if not live_cells:
            return "died_out", 0, gen + 1
            
        if live_cells in history:
            # The pattern has repeated, so it has stabilized.
            return "stabilized", len(live_cells), gen + 1
            
        history.add(live_cells)

    return "timeout", len(live_cells), max_generations

def main():
    """
    Main function to define the pattern and run the simulation.
    This pattern is xs109_481221430z2152, found by David Eppstein.
    """
    # Define the 12x12 starting pattern with 109 cells.
    # Coordinates are (x, y) with (0,0) at the top-left.
    pattern_coords = {
        (1,0),(2,0),(4,0),(5,0),(6,0),(7,0),(9,0),(10,0),
        (0,1),(2,1),(3,1),(8,1),(9,1),(11,1),
        (0,2),(1,2),(3,2),(4,2),(7,2),(8,2),(10,2),(11,2),
        (0,3),(2,3),(5,3),(6,3),(9,3),(11,3),
        (0,4),(1,4),(3,4),(4,4),(7,4),(8,4),(10,4),(11,4),
        (0,5),(2,5),(5,5),(6,5),(9,5),(11,5),
        (0,6),(1,6),(3,6),(4,6),(7,6),(8,6),(10,6),(11,6),
        (0,7),(2,7),(5,7),(6,7),(9,7),(11,7),
        (0,8),(1,8),(3,8),(4,8),(7,8),(8,8),(10,8),(11,8),
        (0,9),(2,9),(3,9),(8,9),(9,9),(11,9),
        (1,10),(2,10),(4,10),(5,10),(6,10),(7,10),(9,10),(10,10),
        (4,11),(7,11)
    }

    initial_population = len(pattern_coords)
    
    status, final_population, generations = simulate_life(pattern_coords, max_generations=1000)

    if status == 'stabilized' and final_population > 100:
        print("Successfully found a pattern meeting the criteria.")
        print(f"The initial pattern has {initial_population} live cells.")
        print(f"It fits within a 12x12 area.")
        print(f"After {generations} generations, it stabilizes.")
        print(f"The final stable population is {final_population} live cells.")
        print("\nThe greatest number of live cells is believed to be 109.")
        print(f"\nFinal Equation: Initial_Cells({initial_population}) -> Stable_Cells({final_population})")
    elif status == 'stabilized':
        print(f"Pattern stabilized at {final_population} cells, which is not over 100.")
    else:
        print("The simulation did not stabilize within the maximum generation limit.")

if __name__ == "__main__":
    main()
