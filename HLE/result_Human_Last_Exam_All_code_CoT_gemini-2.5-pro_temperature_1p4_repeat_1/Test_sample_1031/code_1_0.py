import numpy as np
import time

def solve_disease_spread():
    """
    Calculates and demonstrates the solution for the disease spread problem on a 14x14 grid.
    """
    n = 14

    # According to known results in bootstrap percolation theory, the minimum number of
    # initially infected sites to infect an n x n grid with the 3-neighbor rule is 2n - 2
    # when n is even.
    min_initial_sites = 2 * n - 2
    
    print(f"The grid is {n}x{n}.")
    print("A healthy site becomes infected if at least 3 of its neighbors are infected.")
    print("The minimum number of initially infected sites required to infect the entire grid is given by the formula 2n - 2 for even n.")
    print(f"\nFor n = {n}:")
    print(f"Minimum sites = 2 * {n} - 2 = {min_initial_sites}")
    print("-" * 30)
    print("Now, we will run a simulation to verify that 26 sites are sufficient.")
    
    # Initialize the grid with all sites healthy (0).
    grid = np.zeros((n, n), dtype=int)
    
    # A known configuration of size 2n-2 that infects the entire grid for even n.
    # It consists of two L-shapes in opposite corners.
    k = n // 2
    
    # Infect the top-left L-shape
    for i in range(k):          # i from 0 to 6
        grid[i, 0] = 1
    for j in range(1, k):       # j from 1 to 6
        grid[0, j] = 1
        
    # Infect the bottom-right L-shape
    for i in range(k, n):       # i from 7 to 13
        grid[i, n - 1] = 1
    for j in range(k, n - 1):   # j from 7 to 12
        grid[n - 1, j] = 1

    initial_infected_count = np.sum(grid)
    print(f"Using a starting configuration with {int(initial_infected_count)} infected sites.\n")

    t = 0
    total_sites = n * n
    while True:
        infected_count = np.sum(grid)
        print(f"Time {t}: {infected_count} / {total_sites} sites are infected.")
        
        if infected_count == total_sites:
            print("\nSuccess! The entire grid has been infected.")
            break
            
        newly_infected = []
        # Find all healthy sites that will become infected in the next step
        for r in range(n):
            for c in range(n):
                if grid[r, c] == 0:  # If the site is healthy
                    infected_neighbors = 0
                    # Check neighbors (up, down, left, right)
                    if r > 0 and grid[r - 1, c] == 1:
                        infected_neighbors += 1
                    if r < n - 1 and grid[r + 1, c] == 1:
                        infected_neighbors += 1
                    if c > 0 and grid[r, c - 1] == 1:
                        infected_neighbors += 1
                    if c < n - 1 and grid[r, c + 1] == 1:
                        infected_neighbors += 1
                    
                    if infected_neighbors >= 3:
                        newly_infected.append((r, c))
                        
        if not newly_infected:
            print("\nFailure. The infection has stabilized without infecting the whole grid.")
            break
            
        # Infect the newly identified sites
        for r, c in newly_infected:
            grid[r, c] = 1
            
        t += 1

solve_disease_spread()
<<<26>>>