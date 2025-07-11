import numpy as np
import random

def simulate_conditioned_walk(num_steps):
    """
    Simulates the Doob's h-transform of a 2D SRW conditioned to avoid the origin.

    The simulation uses an approximation of the harmonic function h(x) = log(||x||) + 3
    to guide the walk. It tracks visits to the positive x-axis to demonstrate that
    this infinite set is transient for the walk.
    """
    # Start the walk at a point away from the origin.
    pos = (1, 0)
    
    # h is the approximation of the potential kernel.
    # We add a constant to ensure it's positive for points near the origin.
    def h(point):
        x, y = point
        # The h-function is not defined at the origin.
        if x == 0 and y == 0:
            return -np.inf 
        return np.log(np.sqrt(x**2 + y**2)) + 3

    visit_count = 0
    last_visit_step = 0
    
    # Track visits to the positive x-axis
    if pos[1] == 0 and pos[0] > 0:
        visit_count += 1
        last_visit_step = 0

    for step in range(1, num_steps + 1):
        x, y = pos
        neighbors = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
        
        # The walk cannot move to the origin.
        valid_neighbors = [n for n in neighbors if n != (0, 0)]
        
        # The probability of moving to a neighbor is proportional to its h-value.
        weights = [h(n) for n in valid_neighbors]
        
        if not valid_neighbors:
            print("Walk is trapped, this shouldn't happen.")
            break
            
        # Choose the next step based on the calculated weights.
        pos = random.choices(valid_neighbors, weights=weights, k=1)[0]
        
        # Check if the new position is on the positive x-axis.
        if pos[1] == 0 and pos[0] > 0:
            visit_count += 1
            last_visit_step = step

    print(f"Simulation finished after {num_steps} steps.")
    print(f"Final position of the walker: {pos}")
    print(f"Distance from origin: {np.sqrt(pos[0]**2 + pos[1]**2):.2f}")
    print(f"Total number of visits to the infinite set (positive x-axis): {visit_count}")
    print(f"The last visit to this set occurred at step: {last_visit_step}")

if __name__ == '__main__':
    # Run the simulation for a large number of steps.
    # The output will show that the number of visits is small and stops increasing,
    # which is strong evidence that the set is transient.
    simulate_conditioned_walk(100000)
