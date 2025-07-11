import math
import random

def solve():
    """
    Solves the problem by running a Monte Carlo simulation to find the probability
    that a randomly thrown circle of radius 6 intersects exactly 47 cells of a unit grid.
    """
    R = 6.0
    R2 = R * R
    TARGET_INTERSECTIONS = 47
    
    # A larger number of simulations provides a more accurate result.
    # 200,000 simulations should provide a reasonably stable estimate.
    N_SIMULATIONS = 200000 
    
    # Define the range for cell indices (m, n) to check.
    # A circle with its center in [0,1]x[0,1] and radius 6 can intersect cells
    # with indices m and n from -6 to 6. A slightly larger range [-7, 7] is safe.
    m_range = range(-7, 8)
    n_range = range(-7, 8)
    
    count_target = 0
    
    for _ in range(N_SIMULATIONS):
        # Generate a random center (x, y) in the unit square [0,1]x[0,1]
        x = random.random()
        y = random.random()

        intersections = 0
        for m in m_range:
            for n in n_range:
                # To check if the circumference intersects the cell [m, m+1] x [n, n+1], we check if
                # the cell is neither completely inside the circle's disk nor completely outside.
                
                # 1. Check if the cell is at least partially inside the disk (d_min <= R).
                # d_min is the shortest distance from the center (x,y) to any point in the cell.
                dx_min = 0
                if x < m:
                    dx_min = m - x
                elif x > m + 1:
                    dx_min = x - (m + 1)
                
                dy_min = 0
                if y < n:
                    dy_min = n - y
                elif y > n + 1:
                    dy_min = y - (n + 1)
                
                d_min_sq = dx_min**2 + dy_min**2

                if d_min_sq <= R2:
                    # 2. Check if the cell is not completely contained within the disk (d_max > R).
                    # d_max is the largest distance from the center (x,y) to any point in the cell.
                    dx_max = max(abs(m - x), abs(x - (m + 1)))
                    dy_max = max(abs(n - y), abs(y - (n + 1)))
                    d_max_sq = dx_max**2 + dy_max**2
                    
                    if d_max_sq > R2:
                        intersections += 1
        
        if intersections == TARGET_INTERSECTIONS:
            count_target += 1
            
    probability = count_target / N_SIMULATIONS
    
    print(f"The calculation is based on a Monte Carlo simulation.")
    print(f"Circumference radius R = {R}")
    print(f"Grid mesh size = 1")
    print(f"Target number of intersected cells = {TARGET_INTERSECTIONS}")
    print(f"Total simulation trials N = {N_SIMULATIONS}")
    print(f"Number of trials where intersections = {TARGET_INTERSECTIONS}: {count_target}")
    
    # Outputting the final equation as requested
    print(f"\nFinal Equation for Probability (P):")
    print(f"P = {count_target} / {N_SIMULATIONS}")

    # Outputting the approximate solution with four significant digits
    print(f"\nApproximate probability = {probability:.4g}")
    print(f"<<<{probability:.4g}>>>")

if __name__ == '__main__':
    solve()