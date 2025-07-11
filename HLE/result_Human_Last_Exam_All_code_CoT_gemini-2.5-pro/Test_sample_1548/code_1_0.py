import numpy as np
import math

def simulate_h_transformed_walk(num_steps):
    """
    Simulates a random walk on Z^2 \ {0} with Doob's h-transform,
    using h(x) = log||x|| as the harmonic function.

    The walk is biased to move away from the origin. We test whether the
    x-axis, an infinite set, is transient by counting the number of visits.

    Args:
        num_steps (int): The number of steps to simulate.

    Returns:
        None. Prints the simulation results.
    """
    # Start the walk at (1, 0)
    pos = np.array([1, 0])
    
    # We will check for visits to the x-axis (where y=0)
    # This is our infinite set A. We start on it, so visits start at 1.
    visits_to_x_axis = 1
    
    print(f"Starting simulation for {num_steps} steps.")
    print(f"The chosen infinite set is the x-axis (points where y=0).")
    
    for step in range(num_steps):
        # Define the four neighbors
        neighbors = np.array([
            [pos[0] + 1, pos[1]],
            [pos[0] - 1, pos[1]],
            [pos[0], pos[1] + 1],
            [pos[0], pos[1] - 1]
        ])
        
        weights = []
        valid_neighbors = []
        
        # Calculate weights for each valid neighbor
        for neighbor in neighbors:
            # The walk is on Z^2 \ {0}, so the origin is forbidden.
            if neighbor[0] == 0 and neighbor[1] == 0:
                continue
            
            # The weight is h(neighbor) = log(||neighbor||)
            norm = np.linalg.norm(neighbor)
            weights.append(math.log(norm))
            valid_neighbors.append(neighbor)
            
        # Normalize weights to get probabilities
        total_weight = sum(weights)
        probabilities = [w / total_weight for w in weights]
        
        # Choose the next position based on the probabilities
        next_pos_index = np.random.choice(len(valid_neighbors), p=probabilities)
        pos = valid_neighbors[next_pos_index]
        
        # Check if the new position is on the x-axis
        if pos[1] == 0:
            visits_to_x_axis += 1

    # Calculate the final angle
    final_angle_rad = math.atan2(pos[1], pos[0])
    final_angle_deg = math.degrees(final_angle_rad)

    print("\n--- Simulation Results ---")
    print(f"Total steps taken: {num_steps}")
    # The final equation is not literally an equation, but the result of the count.
    # We output the numbers related to the visit count.
    print(f"Number of visits to the infinite set (x-axis): {visits_to_x_axis}")
    print(f"Final position: ({pos[0]}, {pos[1]})")
    print(f"Final distance from origin: {np.linalg.norm(pos):.2f}")
    print(f"Final angle (degrees): {final_angle_deg:.2f}")

    if visits_to_x_axis < num_steps / 100: # Heuristic check
        print("\nObservation: The number of visits to the x-axis is very small compared to")
        print("the total steps, which supports the conclusion that it is a transient set.")

if __name__ == '__main__':
    # Set the number of steps for the simulation
    # A larger number will show the effect more clearly
    # No "equation" to print numbers from, so we will just run the simulation.
    simulate_h_transformed_walk(20000)
