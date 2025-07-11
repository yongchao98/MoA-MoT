import random
import math

def solve_circumference_problem():
    """
    Calculates the probability that a circumference of a given radius intersects
    a specific number of cells in a unit grid.
    """
    
    # Problem parameters
    R = 6.0
    TARGET_CELL_COUNT = 47
    
    # Simulation parameters
    # A higher number of samples gives a more accurate result but takes longer to run.
    # 1,000,000 samples should provide sufficient accuracy and run in a few minutes.
    NUM_SAMPLES = 1_000_000
    
    R_squared = R * R

    def count_intersected_cells(x, y):
        """
        Counts the number of grid cells intersected by a circle with center (x, y).
        The grid cells are unit squares [n, n+1] x [m, m+1].
        """
        count = 0
        
        # We only need to check cells in a bounding box around the circle.
        # The range is determined by the radius R.
        # n,m from -R-1 to R ensures we cover all possible intersected cells
        # since x,y are in [0,1].
        n_range = range(-int(R) - 1, int(R) + 1)
        m_range = range(-int(R) - 1, int(R) + 1)

        for n in n_range:
            for m in m_range:
                # --- Calculate the squared minimum distance from (x, y) to the cell [n, n+1]x[m, m+1] ---
                # This is the distance to the closest point in the cell.
                dx_min = 0.0
                if x < n:
                    dx_min = n - x
                elif x > n + 1:
                    dx_min = x - (n + 1)
                
                dy_min = 0.0
                if y < m:
                    dy_min = m - y
                elif y > m + 1:
                    dy_min = y - (m + 1)
                
                d_min_sq = dx_min**2 + dy_min**2
                
                # If the closest point is already further than R, the whole cell is, so we can skip.
                if d_min_sq > R_squared:
                    continue
                
                # --- Calculate the squared maximum distance from (x, y) to the cell ---
                # This is the distance to the furthest of the four corners of the cell.
                d_max_sq = max(
                    (n - x)**2 + (m - y)**2,
                    (n + 1 - x)**2 + (m - y)**2,
                    (n - x)**2 + (m + 1 - y)**2,
                    (n + 1 - x)**2 + (m + 1 - y)**2
                )
                
                # The circumference intersects the cell if R is between the min and max distance.
                if d_min_sq <= R_squared <= d_max_sq:
                    count += 1
        return count

    # Run the Monte Carlo simulation
    count_target = 0
    for i in range(NUM_SAMPLES):
        # Generate a random center (x, y) in the unit square [0, 1] x [0, 1]
        rand_x = random.random()
        rand_y = random.random()
        
        num_cells = count_intersected_cells(rand_x, rand_y)
        
        if num_cells == TARGET_CELL_COUNT:
            count_target += 1
            
    # Calculate the final probability
    probability = count_target / NUM_SAMPLES
    
    print(f"A circumference of radius R={R} is thrown at random onto a square grid of mesh 1.")
    print(f"We calculate the probability that it intersects k={TARGET_CELL_COUNT} cells.")
    print("-" * 30)
    print(f"Total number of random throws: {NUM_SAMPLES}")
    print(f"Number of times k={TARGET_CELL_COUNT} intersections occurred: {count_target}")
    print(f"The final equation is: P(N=k) = {count_target} / {NUM_SAMPLES}")
    print("-" * 30)
    # Format the result to four significant digits
    print(f"The approximate probability is: {'%.4g' % probability}")

solve_circumference_problem()